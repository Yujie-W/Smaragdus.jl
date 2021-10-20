module Aeris

using DocStringExtensions: METHODLIST
using PkgUtility: upper_quadratic


# export public functions
export canopy_airspace_fluxes!


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Oct-19: add case of turbulent airspace fluxes using SurfaceFluxes.jl module
# To-do
#     TODO: a method to prescribe airspace fluxes
#
#######################################################################################################################################################################################################
"""
Function to update canopy airspace fluxes. Supported methods are

$(METHODLIST)
"""
function canopy_airspace_fluxes! end


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Oct-19: add case of turbulent airspace fluxes using SurfaceFluxes.jl module
#
#######################################################################################################################################################################################################
"""
"""
canopy_airspace_fluxes!(roughness::FT) where {FT<:AbstractFloat} = (
    return nothing
);


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Oct-20: add function to calculate the canopy length for momentum absorption
# To-do
#     TODO: a function to calculate leaf level drag coefficient?
#
#######################################################################################################################################################################################################
"""
    canopy_length(lai::FT, sai::FT, h::FT, c_d::FT = FT(0.25)) where {FT<:AbstractFloat}

Return the single length scale used to characterize the ability of canopy to absorb momentum, given
- `lai` Leaf area index in `[m² m⁻²]`
- `sai` Stem area index in `[m² m⁻²]`
- `h` Canopy height in `[m]`
- `c_d` Leaf level drag coefficient

# References
- Harman and Finnigan (2007) A simple unified theory for flow in the canopy and roughness sublayer. Boundary-Layer Meteorology 123: 339-363.
- Bonan et al. (2018) Modeling canopy-induced turbulance in the Earth system: A unified parameterization of turbulent exchange within plant canopies and the roughness sublayer (CLM-ml v0).
  Geoscientific Model Development 11: 1467-1496.
"""
function canopy_length(lai::FT, sai::FT, h::FT, c_d::FT = FT(0.25)) where {FT<:AbstractFloat}
    # Harman and Finigan (2007) uses leaf area index for the calculation, and Bonan et al. (2018) extend it to the sum of leaf area index and stem area index.
    return h / (c_d * (lai + sai))
end


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Oct-20: add function to calculate β from Monin-Obukhov length
# To-do
#     TODO: calculate β_n from leaf area index for sparse canopy?
#
#######################################################################################################################################################################################################
"""
    toc_u_ratio(l_c::FT, l_mo::FT, β_n::FT = FT(0.35)) where {FT<:AbstractFloat}

Return the ratio between friction velocity and wind speed at the top of the canopy, given
- `l_c` Single length scale used to characterize the ability of canopy to absorb momentum
- `l_mo` Monin-Obukhov length
- `β_n` Neutral value of the ratio (0.35 by default)

# References
- Bonan et al. (2018) Modeling canopy-induced turbulance in the Earth system: A unified parameterization of turbulent exchange within plant canopies and the roughness sublayer (CLM-ml v0).
  Geoscientific Model Development 11: 1467-1496.
"""
function toc_u_ratio(l_c::FT, l_mo::FT, β_n::FT = FT(0.35)) where {FT<:AbstractFloat}
    @assert l_mo != 0;

    # if Monin-Obukhov length is negative, chose the maximum of the two roots of _β2^2 + 16*l_c/l_mo*β_n^4 * _β2 - β_n^4 = 0
    if l_mo < 0
        _a = FT(1);
        _b = 16 * l_c / l_mo * β_n ^ 4;
        _c = -1 * β_n ^ 4;
        _β2 = upper_quadratic(_a, _b, _c);

        return sqrt(_β2)
    # if l_mo is positive , choose the only real solution of 5*l_c/l_mo * _β^3 + _β - β_n = 0
    else
        _d = 5 * l_c / l_mo;
        _b = β_n;
        _x = (sqrt(3) * sqrt(27*_b^2*_d^4 + 4*_d^3) + 9 * _b *_d^2) ^ (1/3);

        return _x/ (2^(1/3) * 3^(2/3) * _d) - (2/3)^(1/3) / _x
    end
end


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Oct-20: add function to calculate the Schmidt number
#
#######################################################################################################################################################################################################
"""
    schmidt_number(l_c::FT, l_mo::FT) where {FT<:AbstractFloat}

Return Schmidt number, given
- `l_c` Single length scale used to characterize the ability of canopy to absorb momentum
- `l_mo` Monin-Obukhov length

# References
- Harman and Finnigan (2008) Scalar concentration profiles in the canopy and roughness sublayer. Boundary-Layer Meteorology 129: 323-351.
"""
function schmidt_number(l_c::FT, l_mo::FT) where {FT<:AbstractFloat}
    return FT(0.5) + FT(0.3) * tanh(2 * l_c / l_mo)
end


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Oct-20: add function to calculate mixing length of momentum
#
#######################################################################################################################################################################################################
"""
    mixing_length(β::FT, l_c::FT) where {FT<:AbstractFloat}

Return the mixing length of momentum, given
- `β` Ratio between friction velocity and wind speed at the top of the canopy
- `l_c` Single length scale used to characterize the ability of canopy to absorb momentum

# References
- Harman and Finnigan (2007) A simple unified theory for flow in the canopy and roughness sublayer. Boundary-Layer Meteorology 123: 339-363.
"""
function mixing_length(β::FT, l_c::FT) where {FT<:AbstractFloat}
    return 2 * β^3 * l_c
end


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Oct-20: add function to calculate displacement height
#
#######################################################################################################################################################################################################
"""
    displacement_height(β::FT, l_c::FT, h::FT) where {FT<:AbstractFloat}

Return the displacement height within the canopy, given
- `β` Ratio between friction velocity and wind speed at the top of the canopy
- `l_c` Single length scale used to characterize the ability of canopy to absorb momentum
- `h` Canopy height in `[m]`

# References
- Harman and Finnigan (2007) A simple unified theory for flow in the canopy and roughness sublayer. Boundary-Layer Meteorology 123: 339-363.
"""
function displacement_height(β::FT, l_c::FT, h::FT) where {FT<:AbstractFloat}
    return h - β^2 * l_c
end

#=
"""
    wind_speed(wind_toc::FT, u_star_toc::FT, lai::FT, sai::FT, h::FT, z::FT, c_d::FT = FT(0.25)) where {FT<:AbstractFloat}

Return the wind speed within the canopy, given
- `wind_toc` Wind speed at the top of canopy in `[m s⁻¹]`
- `u_star_toc` Friction velocity at the top of canopy in `[m s⁻¹]`
- `lai` Leaf area index in `[m² m⁻²]`
- `sai` Stem area index in `[m² m⁻²]`
- `h` Canopy height in `[m]`
- `z` Height (<= `h`) in `[m]`
- `c_d` Leaf level drag coefficient
"""
function wind_speed(wind_toc::FT, u_star_toc::FT, lai::FT, sai::FT, h::FT, z::FT, c_d::FT = FT(0.25)) where {FT<:AbstractFloat}
    @assert 0 <= z <= h;

    # ratio between friction velocity and wind speed at the top of the canopy
    _β = u_star_toc / wind_toc;

    # calculate mixing length of momentum and displacement height
    _l_m, _ = l_m_and_z_d(_β, lai, sai, h, c_d);

    return wind_toc * exp(_β * (z - h) / _l_m)
end
=#


end # module
