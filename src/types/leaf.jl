"""
$(TYPEDEF)

Struct that contains leaf biophysical traits used to run leaf reflection and transmittance.

# Fields
$(TYPEDFIELDS)
"""
mutable struct LeafBio{FT<:AbstractFloat}
    # parameters that do not change with time
    "Leaf fluorescence quantum efficiency (Fo standard)"
    FQE::FT
    "Dry matter content (dry leaf mass per unit area) `[g cm⁻²]`"
    LMA::FT
    "Leaf mesophyll structural parameter that describes mesophyll reflection and transmittance"
    MESOPHYLL_N::FT
    "Doubling adding layers"
    NDUB::Int

    # variables that change with time
    "Anthocynanin content `[ug cm⁻²]`"
    ant::FT
    "Chlorophyll a and b content `[ug cm⁻²]`"
    cab::FT
    "Carotenoid content `[ug cm⁻²]`"
    car::FT
    "Senescent material fraction `[-]`"
    f_sens::FT
    "Zeaxanthin fraction in Carotenoid (1=all Zeaxanthin, 0=all Violaxanthin) `[-]`"
    f_zeax
    "Relative absorbtion by Chlorophyll `[-]`"
    k_cab::Vector{FT}
    "Relative absorbtion by Chlorophyll+Carotenoid `[-]`"
    k_cabcar::Vector{FT}
    "Equivalent water thickness `[cm]`"
    l_H₂O::FT
    "Fluorescence excitation matrix backwards `[-]`"
    mat_b::Matrix{FT}
    "Fluorescence excitation matrix forwards `[-]`"
    mat_f::Matrix{FT}
    "Broadband thermal reflectance, related to blackbody emittance `[-]`"
    ρ_LW::FT
    "Shortwave leaf reflectance `[-]`"
    ρ_SW::Vector{FT}
    "Broadband thermal transmission, related to blackbody emittance `[-]`"
    τ_LW::FT
    "Shortwave leaf transmission `[-]`"
    τ_SW::Vector{FT}

    # caches to speed up calculations
    "Shortwave absorption, 1 .- ρ_SW .- τ_SW  `[-]`"
    _α_SW::Vector{FT}

    # constructors
    LeafBio{FT}(wls=WaveLengthSet{FT}()) where {FT<:AbstractFloat} = (
        @unpack NΛ, NΛ_SIF, NΛ_SIFE = wls;

        return new{FT}(1, 0.012, 1.4, 10, 0, 40, 10, 0, 0, zeros(FT,NΛ), zeros(FT,NΛ), 0.01, zeros(FT,NΛ_SIF,NΛ_SIFE), zeros(FT,NΛ_SIF,NΛ_SIFE), 0.01,
                       zeros(FT,NΛ), 0.01, zeros(FT,NΛ), zeros(FT,NΛ))
    )
end


"""
$(TYPEDEF)

Struct that contains leaf traits such as plant photosynthesis, hydraulics, and biophysics.

# Fields
$(TYPEDFIELDS)
"""
mutable struct Leaf{FT<:AbstractFloat}
    # parameters that do not change with time

    # variables that change with time
    "Current leaf temperature"
    t::FT

    # caches to speed up calculations
    "Last leaf temperature. If different from t, then make temperature correction"
    _t_old::FT

    # constructors
    Leaf{FT}() where {FT<:AbstractFloat} = new{FT}(298.15, 0.0)
end
