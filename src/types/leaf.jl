#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2021-Nov-04: refactor the structure with constants, variables, and temporary cache
#     2021-Nov-11: add leaf photosynthesis system to the structure
#
#######################################################################################################################################################################################################
mutable struct LeafPhotosynthesis{FT<:AbstractFloat}
    # parameters that do not change with time
    "Fraction of absorbed light used by PSII electron transport reaction `[-]`"
    F_PSII::FT
    "Rate constant for fluorescence `[-]`"
    K_FLUORESCENCE::FT
    "Maximal photochemistry rate constant (all reaction centers open) `[-]`"
    K_PHOTO_MAX::FT
    "Rate constant for thermal dissipation `[-]`"
    K_THERMAL::FT
    "Maximal PSII yield (NPQ=0, all reaction centers open) `[-]`"
    Φ_PSII_MAX::FT

    # prognostic variables that change with time
    "CO₂ partial pressure at leaf chloroplasts, based on Vcmax,c `[Pa]`"
    p_CO₂_chloroplast::FT
    "CO₂ partial pressure at leaf internal airspace, based on Vcmax,i `[Pa]`"
    p_CO₂_internal::FT
    "CO₂ partial pressure at leaf surface, namely external stomatal opening `[Pa]`"
    p_CO₂_surface::FT

    # dignostic variables that change with time
    "Photochemistry rate constant `[-]`"
    k_photo::FT
    "Reversible NPQ rate constant `[-]`"
    k_rev_npq::FT
    "Sustained NPQ rate constant for seasonal changes `[-]`"
    k_sus_npq::FT

    # caches to speed up calculations

    # constructors
    LeafPhotosynthesis{FT}() where {FT<:AbstractFloat} = new{FT}(0.5, 0.05, 4, 0.85, 4/4.9, 20, 30, 40, 4, 0, 0)
end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2021-Aug-04: refactor the structure with BIO_PHYSICS as a field
#     2021-Oct-19: sort variable to prognostic and dignostic catergories
#     2021-Nov-04: add leaf photosynthesis with PHOTOSYNTHESIS as a fields
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
leaf = Leaf{FT}(collect(FT,400:5:2500.1); opti=Smaragdus.OPTI_2021);
leaf = Leaf{FT}(WaveLengthSet{FT}());
leaf = Leaf{FT}(WaveLengthSet{FT}(); opti=Smaragdus.OPTI_2021);
```
"""
mutable struct Leaf{FT<:AbstractFloat}
    # parameters that do not change with time
    "Biophysical parameter structure"
    BIO_PHYSICS::LeafBiophysics{FT}
    "Photosynthetical parameter structure"
    PHOTOSYNTHESIS::LeafPhotosynthesis{FT}

    # prognostic variables that change with time
    "Current leaf temperature"
    t::FT

    # dignostic variables that change with time
    "Saturation H₂O vapor pressure, need to update with temperature and leaf water pressure `[Pa]`"
    p_H₂O_sat::FT

    # caches to speed up calculations
    "Last leaf temperature. If different from t, then make temperature correction"
    _t::FT

    # constructors
    Leaf{FT}(wls::WaveLengthSet{FT}) where {FT<:AbstractFloat} = (
        return new{FT}(LeafBiophysics{FT}(wls), LeafPhotosynthesis{FT}(), 298.15, saturation_vapor_pressure(298.15), 0.0)
    );
end
