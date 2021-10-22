module Folium

using DocStringExtensions: METHODLIST
using LinearAlgebra: mul!
using PkgUtility: numerical∫
using SpecialFunctions: expint
using UnPack: @unpack

# using external types and functions
using ..Emerald: HyperspectralRadiation, Leaf, WaveLengthSet
using ..Radiatio: average_transmittance, energy, photon


# export the public functions
export leaf_spectra!


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Oct-01: rename the function to leaf_spectra! as the function updates not only fluorescence but also reflectance, transmittance, and absorption spectra
#     2021-Oct-22: add another method to prescribe leaf spectra such as transmittance and reflectance from broadband method
#
#######################################################################################################################################################################################################
"""
This function updates leaf level reflectance, transmittance, and fluorescence spectra related parameters. Supported methods are

$(METHODLIST)
"""
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
#     TODO: speed up this function by preallocate memories using a cache structure
#
#######################################################################################################################################################################################################
"""
    leaf_spectra!(leaf::Leaf{FT}, wls::WaveLengthSet{FT}; APAR_car::Bool = true, α::FT=FT(40)) where {FT<:AbstractFloat}

Update leaf reflectance and transmittance spectra, and fluorescence spectrum matrices, given
- `leaf` [`Leaf`](@ref) type struct that contains leaf biophysical parameters
- `wls` [`WaveLengthSet`](@ref) type struct that contain wave length bins
- `APAR_car` If true (default), account carotenoid absorption as APAR; otherwise, APAR is only by chlorophyll
- `α` Optimum angle of incidence (default is 40° as in PROSPECT-D, SCOPE uses 59°)

# Examples
```julia
leaf = Leaf{Float64}();
wls = WaveLengthSet{Float64}();
leaf_spectra!(leaf, wls);
leaf_spectra!(leaf, wls; APAR_car=false);
leaf_spectra!(leaf, wls; APAR_car=false, α=59.0);
```
"""
leaf_spectra!(leaf::Leaf{FT}, wls::WaveLengthSet{FT}; APAR_car::Bool = true, α::FT=FT(40)) where {FT<:AbstractFloat} = (
    BIO = leaf.BIO_PHYSICS;
    @unpack K_ANT, K_BROWN, K_CAB, K_CAR_V, K_CAR_Z, K_CBC, K_H₂O, K_LMA, K_PRO, K_PS, MESOPHYLL_N, NDUB, NR = BIO;
    @unpack IΛ_SIF, IΛ_SIFE, Λ_SIF, Λ_SIFE = wls;

    # calculate the average absorption feature and relative Cab and Car partitions
    BIO._k_all   .= (K_CAB   .* BIO.cab .+                      # chlorophyll absorption
                     K_CAR_V .* BIO.car .* (1 - BIO.f_zeax) .+  # violaxanthin carotenoid absorption
                     K_CAR_Z .* BIO.car .* BIO.f_zeax .+        # zeaxanthin carotenoid absorption
                     K_ANT   .* BIO.ant .+                      # anthocynanin absorption absorption
                     K_BROWN .* BIO.brown .+                    # TODO: needs to be a concentration
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
    BIO.α_SW = 1 .- BIO.τ_SW .- BIO.ρ_SW;

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
    _a   = _x_e .* (1 .+ _y_e .* _y_f) .* _x_f;
    _b   = _x_e .* (_y_e .+ _y_f) .* _x_f;

    BIO.mat_b = _a .* _mat_b + _b .* _mat_f;
    BIO.mat_f = _a .* _mat_f + _b .* _mat_b;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes made to this method
# General
#     2021-Oct-22: add another method to prescribe leaf spectra such as transmittance and reflectance from broadband method
#
#######################################################################################################################################################################################################
"""
    leaf_spectra!(leaf::Leaf{FT}, wls::WaveLengthSet{FT}, ρₚₐᵣ::FT, ρₙᵢᵣ::FT, τₚₐᵣ::FT, τₙᵢᵣ::FT) where {FT<:AbstractFloat}

Update leaf reflectance and transmittance (e.g., prescribe broadband PAR and NIR values), given
- `leaf` [`Leaf`](@ref) type struct that contains leaf biophysical parameters
- `wls` [`WaveLengthSet`](@ref) type struct that contain wave length bins
- `ρₚₐᵣ` Reflectance at PAR region
- `ρₙᵢᵣ` Reflectance at NIR region
- `τₚₐᵣ` Transmittance at PAR region
- `τₙᵢᵣ` Transmittance at NIR region

# Examples
```julia
leaf = Leaf{Float64}();
wls = WaveLengthSet{Float64}();
leaf_spectra!(leaf, wls, 0.1, 0.45, 0.05, 0.25);
```
"""
leaf_spectra!(leaf::Leaf{FT}, wls::WaveLengthSet{FT}, ρₚₐᵣ::FT, ρₙᵢᵣ::FT, τₚₐᵣ::FT, τₙᵢᵣ::FT) where {FT<:AbstractFloat} = (
    BIO = leaf.BIO_PHYSICS;
    @unpack IΛ_NIR, IΛ_PAR = wls;

    BIO.ρ_SW[IΛ_PAR] .= ρₚₐᵣ;
    BIO.ρ_SW[IΛ_NIR] .= ρₙᵢᵣ;
    BIO.τ_SW[IΛ_PAR] .= τₚₐᵣ;
    BIO.τ_SW[IΛ_NIR] .= τₙᵢᵣ;
    BIO.α_SW = 1 .- BIO.τ_SW .- BIO.ρ_SW;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Oct-22: add function to compute leaf level PAR and APAR
# To do
#     TODO: use cache to speed this up
#
#######################################################################################################################################################################################################
"""
    leaf_PAR(leaf::Leaf{FT}, wls::WaveLengthSet{FT}, rad::HyperspectralRadiation{FT}; APAR_car::Bool = true) where {FT<:AbstractFloat}

Return leaf level PAR and APAR, given
- `leaf` [`Leaf`](@ref) type struct that contains leaf biophysical parameters
- `wls` [`WaveLengthSet`](@ref) type struct that contains wave length bins
- `rad` [`HyperspectralRadiation`](@ref) type struct that contains incoming radiation information
- `APAR_car` If true (default), account carotenoid absorption as APAR; otherwise, APAR is only by chlorophyll

---
# Examples
```julia
leaf = Leaf{Float64}();
wls = WaveLengthSet{Float64}();
rad = HyperspectralRadiation{Float64}();
par,apar = leaf_PAR(leaf, wls, rad);
par,apar = leaf_PAR(leaf, wls, rad; APAR_car=false);
```
"""
function leaf_PAR(leaf::Leaf{FT}, wls::WaveLengthSet{FT}, rad::HyperspectralRadiation{FT}; APAR_car::Bool = true) where {FT<:AbstractFloat}
    BIO = leaf.BIO_PHYSICS;
    @unpack α_cab, α_cabcar, α_SW = BIO;
    @unpack e_direct, e_diffuse = rad;
    @unpack IΛ_PAR, ΔΛ_PAR, Λ_PAR = wls;

    # APAR absorption feature
    _α = (APAR_car ? view(α_cabcar, IΛ_PAR) : view(α_cab, IΛ_PAR));

    # PAR energy from direct  and diffuse light
    _e_par_dir  = view(e_direct , IΛ_PAR) .* view(α_SW, IΛ_PAR);
    _e_par_diff = view(e_diffuse, IΛ_PAR) .* view(α_SW, IΛ_PAR);
    _par_dir  = photon.(Λ_PAR, _e_par_dir );
    _par_diff = photon.(Λ_PAR, _e_par_diff);

    # absorbed PAR energy from direct and diffuse light
    _apar_dir  = _α .* _par_dir;
    _apar_diff = _α .* _par_diff;

    # total PAR and APAR in μmol photons m⁻² s⁻¹
    _∑par_dir   = numerical∫(_par_dir  , ΔΛ_PAR);
    _∑par_diff  = numerical∫(_par_diff , ΔΛ_PAR);
    _∑apar_dir  = numerical∫(_apar_dir , ΔΛ_PAR);
    _∑apar_diff = numerical∫(_apar_diff, ΔΛ_PAR);

    return 1000 * (_∑par_dir + _∑par_diff), 1000 * (_∑apar_dir + _∑apar_diff)
end


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Jul-08: add leaf level SIF simulation
#     2021-Jul-08: use mat_b and mat_f for SIF at backward and forward directions
#     2021-Aug-05: add option to sumulate SIF in photon to photon mode
#     2021-Oct-22: refactor the function to leaf_SIF to return the SIFs directly
#
#######################################################################################################################################################################################################
"""
    leaf_SIF(leaf::Leaf{FT}, wls::WaveLengthSet{FT}, rad::HyperspectralRadiation{FT}, ϕ::FT = FT(0.01); ϕ_photon::Bool = true) where {FT<:AbstractFloat}

Return the leaf level SIF at backward and forward directions, given
- `leaf` [`Leaf`](@ref) type struct that contains leaf biophysical parameters
- `wls` [`WaveLengthSet`](@ref) type struct that contains wave length bins
- `rad` [`HyperspectralRadiation`](@ref) type struct that contains incoming radiation information
- `ϕ` Fluorescence quantum yield
- `ϕ_photon` If true (default), convert photon to photon when computing SIF; otherwise, convert energy to energy

---
# Examples
```julia
leaf = Leaf{Float64}();
wls = WaveLengthSet{Float64}();
rad = HyperspectralRadiation{Float64}();
sif_b,sif_f = leaf_SIF(leaf, wls, rad, 0.01);
sif_b,sif_f = leaf_SIF(leaf, wls, rad, 0.01; ϕ_photon=false);
```
"""
function leaf_SIF(leaf::Leaf{FT}, wls::WaveLengthSet{FT}, rad::HyperspectralRadiation{FT}, ϕ::FT = FT(0.01); ϕ_photon::Bool = true) where {FT<:AbstractFloat}
    BIO = leaf.BIO_PHYSICS;
    @unpack mat_b, mat_f = BIO;
    @unpack e_direct, e_diffuse = rad;
    @unpack IΛ_SIFE, ΔΛ_SIFE, Λ_SIF, Λ_SIFE = wls;

    # calculate the excitation energy and photons
    _e_excitation = (view(e_direct , IΛ_SIFE) .+ view(e_diffuse, IΛ_SIFE)) .* ΔΛ_SIFE;

    # convert energy to energy using the matrices
    if !ϕ_photon
        _sif_b = mat_b * _e_excitation * ϕ / FT(pi);
        _sif_f = mat_f * _e_excitation * ϕ / FT(pi);

        return _sif_b, _sif_f
    end;

    # convert energy to photon
    _phot_excitation = photon.(Λ_SIFE, _e_excitation);

    # convert photon to photon using the matrices
    _phot_b = mat_b * _phot_excitation * ϕ / FT(pi);
    _phot_f = mat_f * _phot_excitation * ϕ / FT(pi);

    # convert photon to back to energy
    _sif_b = energy.(Λ_SIFE, _phot_b);
    _sif_f = energy.(Λ_SIFE, _phot_f);

    return _sif_b, _sif_f
end


end # module
