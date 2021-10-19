module Aeris

using DocStringExtensions: METHODLIST
using SurfaceFluxes: DGScheme, surface_conditions


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


end # module
