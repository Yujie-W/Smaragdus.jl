#=
#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2021-Sep-30: define this structure with no constructors
#     2021-Oct-19: add soil and air variables
#     2021-Oct-19: sort variable to prognostic and dignostic catergories
#
#######################################################################################################################################################################################################
"""
$(TYPEDEF)

Struct that contains environmental conditions, such as soil moisture and atmospheric vapor pressure. Note that this structure is designed to be containers to interact with other CliMA modules and to
    prescribe values.

# Fields
$(TYPEDFIELDS)

---
# Examples
```julia
;
```
"""
mutable struct SoilAir{FT<:AbstractFloat}
    # parameters that do not change with time
    "Total area of the soil/air interface `[m²]`"
    AREA::FT
    "Soil color class used for soil albedo calculations"
    COLOR::Int
    "Number of air layers"
    N_AIR::Int
    "Number of soil layers"
    N_SOIL::Int
    "Soil moisture retention curve"
    VC::Union{Vector{BrooksCorey{FT}}, Vector{VanGenuchten{FT}}}
    "Z profile for air `[m]`"
    Z_AIR::Vector{FT}
    "Z profile for soil `[m]`"
    Z_SOIL::Vector{FT}
    "ΔZ profile for air `[m]`"
    ΔZ_AIR::Vector{FT}
    "ΔZ profile for soil `[m]`"
    ΔZ_SOIL::Vector{FT}

    # prognostic variables that change with time
    "CO₂ partial pressure at different air layers `[Pa]`"
    p_CO₂::Vector{FT}
    "H₂O partial pressure at different air layers `[Pa]`"
    p_H₂O::Vector{FT}
    "Temperature at different air layers `[K]`"
    t_air::Vector{FT}
    "Temperature at different soil layers `[K]`"
    t_soil::Vector{FT}
    "Wind speed (total) `[m s⁻¹]`"
    wind::Vector{FT}
    "Wind speed (vertical) `[m s⁻¹]`"
    wind_z::Vector{FT}
    "Soil water content at different soil layers `[m³ m⁻³]`"
    θ::Vector{FT}

    # dignostic variables that change with time
    "Saturated vapor pressure at different air layers `[Pa]`"
    p_H₂O_sat::Vector{FT}
    "Relative humidity at different air layers `[-]`"
    rh::Vector{FT}
    "Vapor pressure deficit at different air layers `[Pa]`"
    vpd::Vector{FT}
    "Soil water potential at different soil layers `[MPa]`"
    ψ::Vector{FT}

    # caches to speed up calculations

    # constructors
    SoilAir{FT}(z_soil::Vector{FT}, z_air::Vector{FT}, area::FT=FT(1)) where {FT<:AbstractFloat} = (
        _n_air  = length(z_air) - 1;
        _n_soil = length(z_soil) - 1;
        _p_sats = saturation_vapor_pressure(T_25()) * ones(_n_air);
        _vcs    = [VanGenuchten{FT}("Silt") for _i in 1:_n_soil];

        return new{FT}(area, 1, _n_air, _n_soil, _vcs, z_air, z_soil, diff(z_air), diff(z_soil), 41 * ones(_n_air),  0.5 * _p_sats, T_25() * ones(_n_air), T_25() * ones(_n_soil), 2 * ones(_n_air),
                       1 * ones(_n_air), [_vc.Θ_SAT for _vc in _vcs], _p_sats, 0.5 * ones(_n_air), 0.5 * _p_sats, zeros(_n_soil))
    );
end
=#
