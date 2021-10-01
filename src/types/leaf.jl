#######################################################################################################################################################################################################
# Changes to this structure
# General
#     2021-Aug-04: refactor the structure with constants, variables, and temporary cache
#     2021-Aug-04: add concentrations and characteristic curves altogether
#     2021-Aug-10: add CBC and PRO supoort
#     2021-Agu-10: add constructors within the structure rather than initialize it externally
#     2021-Sep-30: rename LeafBio to LeafBiophysics to be more specific
#
#######################################################################################################################################################################################################
"""
$(TYPEDEF)

Struct that contains leaf biophysical traits used to run leaf reflection and transmittance.

# Fields
$(TYPEDFIELDS)

---
# Examples
```julia
lbio = Emerald.LeafBiophysics{FT}();
lbio = Emerald.LeafBiophysics{FT}(collect(FT,400:5:2500.1));
lbio = Emerald.LeafBiophysics{FT}(collect(FT,400:5:2500.1); opti=Emerald.OPTI_2021);
lbio = Emerald.LeafBiophysics{FT}(WaveLengthSet{FT}());
lbio = Emerald.LeafBiophysics{FT}(WaveLengthSet{FT}(); opti=Emerald.OPTI_2021);
```
"""
mutable struct LeafBiophysics{FT<:AbstractFloat}
    # parameters that do not change with time
    "Leaf fluorescence quantum efficiency (Fo standard)"
    FQE::FT
    "Dry matter content (dry leaf mass per unit area) `[g cm⁻²]`"
    LMA::FT
    "Carbon-based constituents in LMA `[g cm⁻²]`"
    CBC::FT
    "Protein content in LMA (PRO = LMA - CBC) `[g cm⁻²]`"
    PRO::FT
    "Specific absorption coefficients of anthocynanin `[-]`"
    K_ANT::Vector{FT}
    "Specific absorption coefficients of chlorophyll a and b `[-]`"
    K_CAB::Vector{FT}
    "Specific absorption coefficients of violaxanthin carotenoid `[-]`"
    K_CAR_V::Vector{FT}
    "Specific absorption coefficients of zeaxanthin carotenoid `[-]`"
    K_CAR_Z::Vector{FT}
    "Specific absorption coefficients of carbon-based constituents `[-]`"
    K_CBC::Vector{FT}
    "Specific absorption coefficients of water `[-]`"
    K_H₂O::Vector{FT}
    "Specific absorption coefficients of dry matter `[-]`"
    K_LMA::Vector{FT}
    "Specific absorption coefficients of protein `[-]`"
    K_PRO::Vector{FT}
    "Specific absorption coefficients of PS I and II `[-]`"
    K_PS::Vector{FT}
    "Specific absorption coefficients of senescent material `[-]`"
    K_SENES::Vector{FT}
    "Leaf mesophyll structural parameter that describes mesophyll reflection and transmittance"
    MESOPHYLL_N::FT
    "Doubling adding layers"
    NDUB::Int
    "Refractive index `[-]`"
    NR::Vector{FT}

    # variables that change with time
    "Anthocynanin content `[ug cm⁻²]`"
    ant::FT
    "Chlorophyll a and b content `[ug cm⁻²]`"
    cab::FT
    "Carotenoid content `[ug cm⁻²]`"
    car::FT
    "Senescent material (brown pigments) fraction `[-]`"
    f_senes::FT
    "Zeaxanthin fraction in Carotenoid (1=all Zeaxanthin, 0=all Violaxanthin) `[-]`"
    f_zeax
    "Equivalent water thickness `[cm]`"
    l_H₂O::FT
    "Fluorescence excitation matrix backwards `[-]`"
    mat_b::Matrix{FT}
    "Fluorescence excitation matrix forwards `[-]`"
    mat_f::Matrix{FT}
    "Relative absorption by Chlorophyll `[-]`"
    α_cab::Vector{FT}
    "Relative absorption by Chlorophyll+Carotenoid `[-]`"
    α_cabcar::Vector{FT}
    "Broadband thermal reflectance, related to blackbody emittance `[-]`"
    ρ_LW::FT
    "Shortwave leaf reflectance `[-]`"
    ρ_SW::Vector{FT}
    "Broadband thermal transmission, related to blackbody emittance `[-]`"
    τ_LW::FT
    "Shortwave leaf transmission `[-]`"
    τ_SW::Vector{FT}

    # caches to speed up calculations
    "Specific absorption coefficients of all materials"
    _k_all::Vector{FT}
    "Shortwave absorption, 1 .- ρ_SW .- τ_SW  `[-]`"
    _α_SW::Vector{FT}

    # constructors
    LeafBiophysics{FT}(swl::Vector{FT}=FT.(WAVELENGTHS); opti::String=OPTI_2021) where {FT<:AbstractFloat} = LeafBiophysics{FT}(WaveLengthSet{FT}(swl; opti=opti); opti=opti)
    LeafBiophysics{FT}(wls::WaveLengthSet{FT}; opti::String=OPTI_2021) where {FT<:AbstractFloat} = (
        @unpack NΛ, NΛ_SIF, NΛ_SIFE, SΛ = wls;

        # read data from the MAT file
        _opti    = matread(opti)["optipar"];
        __nr     = _opti["nr"  ];
        __Klma   = _opti["Kdm" ];
        __Kcab   = _opti["Kab" ];
        __Kant   = _opti["Kant"];
        __Kh2o   = _opti["Kw"  ];
        __Ksenes = _opti["Ks"  ];
        __Kps    = _opti["phi" ];
        __KcaV   = _opti["KcaV"];
        __KcaZ   = _opti["KcaZ"];
        __lambda = _opti["wl"  ];

        # read protein and carbon-based constituents if exist (not in file OPTI_2017)
        if opti == OPTI_2017
            __Kpro = __Klma;
            __Kcbc = __Klma;
        else
            __Kpro = _opti["Kp"  ];
            __Kcbc = _opti["Kcbc"];
        end;

        # create data to parse
        _nr     = zeros(FT, NΛ);
        _Klma   = zeros(FT, NΛ);
        _Kcab   = zeros(FT, NΛ);
        _Kant   = zeros(FT, NΛ);
        _Kh2o   = zeros(FT, NΛ);
        _Ksenes = zeros(FT, NΛ);
        _Kps    = zeros(FT, NΛ);
        _KcaV   = zeros(FT, NΛ);
        _KcaZ   = zeros(FT, NΛ);
        _Kpro   = zeros(FT, NΛ);
        _Kcbc   = zeros(FT, NΛ);

        # fill in the data arrays
        @inbounds for _i in 1:NΛ
            _wo = findall( (__lambda.>=SΛ[_i]) .& (__lambda.<SΛ[_i+1]) );
            if length(_wo) < 1
                @warn "Warning, some wavelengths out of bounds $(string(SΛ[_i]))";
            end

            _nr[_i]     = nanmean(__nr[_wo]    );
            _Klma[_i]   = nanmean(__Klma[_wo]  );
            _Kcab[_i]   = nanmean(__Kcab[_wo]  );
            _Kant[_i]   = nanmean(__Kant[_wo]  );
            _Kh2o[_i]   = nanmean(__Kh2o[_wo]  );
            _Ksenes[_i] = nanmean(__Ksenes[_wo]);
            _Kps[_i]    = nanmean(__Kps[_wo]   );
            _KcaV[_i]   = nanmean(__KcaV[_wo]  );
            _KcaZ[_i]   = nanmean(__KcaZ[_wo]  );
            _Kpro[_i]   = nanmean(__Kpro[_wo]  );
            _Kcbc[_i]   = nanmean(__Kcbc[_wo]  );
        end;

        return new{FT}(1, 0.012, 0, 0, _Kant, _Kcab, _KcaV, _KcaZ, _Kcbc, _Kh2o, _Klma, _Kpro, _Kps, _Ksenes, 1.4, 10, _nr, 0, 40, 10, 0, 0, 0.01, zeros(FT,NΛ_SIF,NΛ_SIFE), zeros(FT,NΛ_SIF,NΛ_SIFE),
                       zeros(FT,NΛ), zeros(FT,NΛ), 0.01, zeros(FT,NΛ), 0.01, zeros(FT,NΛ), zeros(FT,NΛ), zeros(FT,NΛ))
    )
end


#######################################################################################################################################################################################################
# Changes to this structure
# General
#     2021-Aug-04: refactor the structure with BIO_PHYSICS as a field
# To do
#     TODO: add leaf physiological parameters as a field well
#     TODO: add leaf hydraulics as a field as well
#     TODO: link leaf water content to BIO_PHYSICS.l_H₂O
#
#######################################################################################################################################################################################################
"""
$(TYPEDEF)

Struct that contains leaf traits such as plant photosynthesis, hydraulics, and biophysics.

# Fields
$(TYPEDFIELDS)

---
# Examples
```julia
leaf = Leaf{FT}();
leaf = Leaf{FT}(collect(FT,400:5:2500.1));
leaf = Leaf{FT}(collect(FT,400:5:2500.1); opti=Emerald.OPTI_2021);
leaf = Leaf{FT}(WaveLengthSet{FT}());
leaf = Leaf{FT}(WaveLengthSet{FT}(); opti=Emerald.OPTI_2021);
```
"""
mutable struct Leaf{FT<:AbstractFloat}
    # parameters that do not change with time
    "Biophysical parameter structure"
    BIO_PHYSICS::LeafBiophysics{FT}

    # variables that change with time
    "Current leaf temperature"
    t::FT

    # caches to speed up calculations
    "Last leaf temperature. If different from t, then make temperature correction"
    _t::FT

    # constructors
    Leaf{FT}(wls::WaveLengthSet{FT}; opti::String=OPTI_2021) where {FT<:AbstractFloat} = new{FT}(LeafBiophysics{FT}(wls; opti=opti), 298.15, 0.0)
    Leaf{FT}(swl::Vector{FT}=FT.(WAVELENGTHS); opti::String=OPTI_2021) where {FT<:AbstractFloat} = new{FT}(LeafBiophysics{FT}(swl; opti=opti), 298.15, 0.0)
end
