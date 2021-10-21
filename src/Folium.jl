module Folium

using SpecialFunctions: expint
using UnPack: @unpack

# using external types and functions
using ..Emerald: Leaf, WaveLengthSet
using ..Radiatio: average_transmittance


# export the public functions
export leaf_spectra!


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Oct-01: rename the function to leaf_spectra! as the function updates not only fluorescence but also reflectance, transmittance, and absorption spectra
# Bug fix
# To do
#     TODO: add another method to prescribe leaf spectra such as transmittance and reflectance from broadband method
#
#######################################################################################################################################################################################################
function leaf_spectra! end


#######################################################################################################################################################################################################
#
# Changes made to this method
# General
#     2020-Mar-30: account for carotenoid absorption as APAR as well as chlorophyll
#     2020-Mar-31: use 40° rather than 59° for _τ_α calculation (following PROSPECT-D)
#     2021-Aug-07: replace function `expint` with that from SpecialFunctions
#     2021-Oct-21: add α to input parameters so that one can roll back to 50° for _τ_α calculation
# Bug fix
#     2021-Aug-06: If BIO.CBC and BIO.PRO are not zero, they are accounted for twice in BIO.LMA, thus the spectrum from LMA need to subtract the contribution from CBC and PRO
# To do
#     TODO: make brown pigment and absoption curve more general using realistic units
#     TODO: add References for this methods
#
#######################################################################################################################################################################################################
leaf_spectra!(leaf::Leaf{FT}, wls::WaveLengthSet{FT}; APAR_car::Bool = true, α::FT=FT(40)) where {FT<:AbstractFloat} = (
    BIO = leaf.BIO_PHYSICS;
    @unpack K_ANT, K_BROWN, K_CAB, K_CAR_V, K_CAR_Z, K_CBC, K_H₂O, K_LMA, K_PRO, K_PS, MESOPHYLL_N, NDUB, NR = BIO;
    @unpack IΛ_SIF, IΛ_SIFE, Λ_SIF, Λ_SIFE = wls;

    # calculate the average absorption feature and relative Cab and Car partitions
    BIO._k_all   .= (K_CAB   .* BIO.cab .+                      # chlorophyll absorption
                     K_CAR_V .* BIO.car .* (1 - BIO.f_zeax) .+  # violaxanthin carotenoid absorption
                     K_CAR_Z .* BIO.car .* BIO.f_zeax .+        # zeaxanthin carotenoid absorption
                     K_ANT   .* BIO.ant .+                      # anthocynanin absorption absorption
                     K_BROWN .* BIO.brown .+                    # TODO: cannot be a fraction, needs to be a concentration
                     K_H₂O   .* BIO.l_H₂O .+                    # water absorption
                     K_CBC   .* BIO.CBC .+                      # carbon-based constituents absorption
                     K_PRO   .* BIO.PRO .+                      # protein absorption
                     K_LMA   .* (BIO.LMA - BIO.CBC - BIO.PRO)   # dry mass absorption (if some remained)
                    ) ./ MESOPHYLL_N;
    BIO.α_cab    .= (K_CAB .* BIO.cab) ./ BIO._k_all ./ MESOPHYLL_N;
    BIO.α_cabcar .= (K_CAB .* BIO.cab .+ K_CAR_V .* BIO.car .* (1 - BIO.f_zeax) .+ K_CAR_Z .* BIO.car .* BIO.f_zeax) ./ BIO._k_all ./ MESOPHYLL_N;

    # calculate the reflectance and transmittance at the interfaces of one layer
    _τ   = (1 .- BIO._k_all) .* exp.(-BIO._k_all) .+ BIO._k_all .^ 2 .* expint.(BIO._k_all .+ eps(FT));
    _τ_α = average_transmittance.(α, NR);
    _ρ_α = 1 .- _τ_α;
    _τ₁₂ = average_transmittance.(FT(90), NR);
    _ρ₁₂ = 1 .- _τ₁₂;
    _τ₂₁ = _τ₁₂ ./ (NR .^ 2);
    _ρ₂₁ = 1 .- _τ₂₁;

    # top surface side
    _denom = 1 .- (_τ .* _ρ₂₁) .^ 2;
    _τ_top = _τ_α .* _τ .* _τ₂₁ ./ _denom;
    _ρ_top = _ρ_α .+ _ρ₂₁ .* _τ .* _τ_top;

    # bottom surface side
    _τ_bottom = _τ₁₂ .* _τ .* _τ₂₁ ./ _denom;
    _ρ_bottom = _ρ₁₂ .+ _ρ₂₁ .* _τ .* _τ_bottom;

    # calculate the reflectance and transmittance at the interfaces of N layer
    _d     = sqrt.((1 .+ _ρ_bottom .+ _τ_bottom) .* (1 .+ _ρ_bottom .- _τ_bottom) .* (1 .- _ρ_bottom .+ _τ_bottom) .* (1 .- _ρ_bottom .- _τ_bottom));
    _ρ²    = _ρ_bottom .^ 2;
    _τ²    = _τ_bottom .^ 2;
    _a     = (1 .+ _ρ² .- _τ² .+ _d) ./ (2 .* _ρ_bottom);
    _b     = (1 .- _ρ² .+ _τ² .+ _d) ./ (2 .* _τ_bottom);
    _bⁿ⁻¹  = _b .^ (MESOPHYLL_N - 1);
    _b²ⁿ⁻² = _bⁿ⁻¹ .^ 2;
    _a²    = _a .^ 2;
    _denom = _a² .* _b²ⁿ⁻² .- 1;
    _ρ_sub = _a .* (_b²ⁿ⁻² .- 1) ./ _denom;
    _τ_sub = _bⁿ⁻¹ .* (_a² .- 1) ./ _denom;

    # avoid case of zero absorption
    _j = findall(_ρ_bottom .+ _τ_bottom .>= 1);
    _τ_sub[_j] = _τ_bottom[_j] ./ (_τ_bottom[_j] + (1 .- _τ_bottom[_j]) * (MESOPHYLL_N - 1));
    _ρ_sub[_j] = 1 .- _τ_sub[_j];

    # reflectance & transmittance of the leaf: combine top layer with next N-1 layers
    _denom    = 1 .- _ρ_sub .* _ρ_bottom;
    BIO.τ_SW  = _τ_top .* _τ_sub ./ _denom;
    BIO.ρ_SW  = _ρ_top .+ _τ_top .* _ρ_sub .* _τ_bottom ./ _denom;
    BIO._α_SW = 1 .- BIO.τ_SW .- BIO.ρ_SW;

    # Doubling method used to calculate fluoresence is now only applied to the part of the leaf where absorption takes place, that is, the part exclusive of the leaf-air interfaces.
    # The reflectance (rho) and transmittance (tau) of this part of the leaf are now determined by "subtracting" the interfaces.
    # CF Note: All of the below takes about 10 times more time than the RT above. Need to rething speed and accuracy. (10nm is bringing it down a lot!)
    _ρ_b = (BIO.ρ_SW .- _ρ_α) ./ (_τ_α .* _τ₂₁ .+ (BIO.ρ_SW - _ρ_α) .* _ρ₂₁);
    _tt1 = _τ_α .* _τ₂₁;
    _tt2 = BIO.τ_SW .* (1 .- _ρ_b .* _ρ₂₁);
    _z   = _tt2 ./ _tt1;
    _tt1 = _ρ_b - _ρ₂₁ .* _z .^ 2;
    _tt2 = 1 .- (_ρ₂₁.* _z) .^ 2;
    _ρ   = max.(0, _tt1 ./ _tt2);
    _tt1 = 1 .- _ρ_b .* _ρ₂₁;
    _τ   = _tt1 ./ _tt2 .* _z;

    # Derive Kubelka-Munk s and k
    _i      = findall((_ρ .+ _τ) .< 1);
    _j      = findall((_ρ .+ _τ) .> 1);
    _d[_i] .= sqrt.((1 .+ _ρ[_i] .+ _τ[_i]) .* (1 .+ _ρ[_i] .- _τ[_i]) .* (1 .- _ρ[_i] .+ _τ[_i]) .*  (1 .- _ρ[_i] .- _τ[_i]));
    _a[_i] .= (1 .+ _ρ[_i] .^ 2 .- _τ[_i] .^ 2 .+ _d[_i]) ./ (2 .* _ρ[_i]);
    _b[_i] .= (1 .- _ρ[_i] .^ 2 .+ _τ[_i] .^ 2 .+ _d[_i]) ./ (2 .* _τ[_i]);
    _a[_j] .= 1;
    _b[_j] .= 1;

    _i      = findall((_a .> 1) .& (_a .!= Inf));
    _s      = _ρ ./ _τ;
    _k      = log.(_b);
    _s[_i] .= 2 .* _a[_i] ./ (_a[_i] .^ 2 .- 1) .* log.(_b[_i]);
    _k[_i] .= (_a[_i] .- 1) ./ (_a[_i] .+ 1) .* log.(_b[_i]);
    _k_chl  = (APAR_car ? BIO.α_cabcar : BIO.α_cab) .* _k;

    # indices of WLE and WLF within wlp
    _ϵ       = FT(2) ^ -NDUB;
    _τ_e     = 1 .- (_k[IΛ_SIFE] .+ _s[IΛ_SIFE]) * _ϵ;
    _τ_f     = 1 .- (_k[IΛ_SIF] .+ _s[IΛ_SIF]) * _ϵ;
    _ρ_e     = _s[IΛ_SIFE] * _ϵ;
    _ρ_f     = _s[IΛ_SIF] * _ϵ;
    _sigmoid = 1 ./ (1 .+ exp.(-Λ_SIF ./ 10) .* exp.(Λ_SIFE' ./ 10));
    _mat_f   = K_PS[IΛ_SIF] .* _ϵ ./ 2 .* _k_chl[IΛ_SIFE]' .* _sigmoid;
    _mat_b   = K_PS[IΛ_SIF] .* _ϵ ./ 2 .* _k_chl[IΛ_SIFE]' .* _sigmoid;

    # Doubling adding routine
    _1_h = ones(FT, 1, length(_τ_e));
    _1_v = ones(FT, length(_τ_f), 1);
    for i in 1:NDUB
        _x_e     = _τ_e ./ (1 .- _ρ_e .^ 2);
        _x_f     = _τ_f ./ (1 .- _ρ_f .^ 2);
        _τ_e_n   = _τ_e .* _x_e;
        _τ_f_n   = _τ_f .* _x_f;
        _ρ_e_n   = _ρ_e .* (1 .+ _τ_e_n);
        _ρ_f_n   = _ρ_f .* (1 .+ _τ_f_n);
        _a₁₁     = _x_f * _1_h .+ _1_v * _x_e';
        _a₁₂     = (_x_f * _x_e') .* (_ρ_f * _1_h .+ _1_v * _ρ_e');
        _a₂₁     = 1 .+ (_x_f * _x_e') .* (1 .+ _ρ_f * _ρ_e');
        _a₂₂     = (_x_f .* _ρ_f) * _1_h .+ _1_v * (_x_e.*_ρ_e)';
        _mat_f_n = _mat_f .* _a₁₁ .+ _mat_b .* _a₁₂;
        _mat_b_n = _mat_b .* _a₂₁ .+ _mat_f .* _a₂₂;
        _τ_e     = _τ_e_n;
        _ρ_e     = _ρ_e_n;
        _τ_f     = _τ_f_n;
        _ρ_f     = _ρ_f_n;
        _mat_f   = _mat_f_n;
        _mat_b   = _mat_b_n;
    end;

    # This reduced red SIF quite a bit in backscatter, not sure why.
    _ρ_b = _ρ .+ _τ .^ 2 .* _ρ₂₁ ./ (1 .- _ρ .* _ρ₂₁);
    _x_e = _1_v * (_τ_α[IΛ_SIFE] ./ (1 .- _ρ₂₁[IΛ_SIFE] .* _ρ_b[IΛ_SIFE]))';
    _x_f = _τ₂₁[IΛ_SIF] ./ (1 .- _ρ₂₁[IΛ_SIF] .* _ρ_b[IΛ_SIF]) * _1_h;
    _y_e = _1_v * (_τ[IΛ_SIFE] .* _ρ₂₁[IΛ_SIFE] ./ (1 .- _ρ[IΛ_SIFE] .* _ρ₂₁[IΛ_SIFE]))';
    _y_f = _τ[IΛ_SIF] .* _ρ₂₁[IΛ_SIF] ./ (1 .- _ρ[IΛ_SIF] .* _ρ₂₁[IΛ_SIF]) * _1_h;
    _a   = _x_e .* (1 .+ _y_e.*_y_f) .* _x_f;
    _b   = _x_e .* (_y_e .+ _y_f) .* _x_f;

    BIO.mat_b = _a .* _mat_b + _b .* _mat_f;
    BIO.mat_f = _a .* _mat_f + _b .* _mat_b;

    return nothing
);


end
