module Soli

using ConstrainedRootSolvers: ReduceStepMethodND, SolutionToleranceND, find_peak
using DocStringExtensions: METHODLIST
using ..Emerald: VanGenuchten
using UnPack: @unpack

import ..Emerald: BrooksCorey


# export public functions
export soil_ψ_25


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Sep-30: move this function out of BrooksCorey struct as an external method for the constructor
#
#######################################################################################################################################################################################################
"""
    BrooksCorey{FT}(name::String)

A constructor for BrooksCorey to create BrooksCorey type soil from VanGenuchten type
"""
BrooksCorey{FT}(vg::VanGenuchten{FT}) where {FT<:AbstractFloat} = (
    # generate data to fit
    _Θs   = range(vg.Θ_RES+FT(1e-2); stop=vg.Θ_SAT-FT(1e-2), length=30);
    _Ψ_vG = -1 .* soil_ψ_25.([vg], _Θs);

    # function to fit BrooksCorey parameters
    @inline _fit(x) = (
        _bc   = BrooksCorey{FT}(vg.TYPE, x[1], x[2], vg.Θ_SAT, vg.Θ_RES);
        _Ψ_BC = -1 .* soil_ψ_25.([_bc], _Θs);
        _diff = sum( (log.(_Ψ_BC) .- log.(_Ψ_vG)) .^ 2 );
        return -_diff
    );

    _st  = SolutionToleranceND{FT}([1e-3, 1e-6], 30);
    _ms  = ReduceStepMethodND{FT}(x_mins = FT[1e-3, 1e-6], x_maxs = FT[ 100, 1000], x_inis = [(2*vg.N-1) / (vg.N-1), 1 / (vg.Α)], Δ_inis = FT[0.1, 1e-3]);
    _sol = find_peak(_fit, _ms, _st);

    return BrooksCorey{FT}(vg.TYPE, _sol[1], _sol[2], vg.Θ_SAT, vg.Θ_RES);
)


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Sep-30: create this function to work with two soil types using either VanGenuchten or BrooksCorey function
#
#######################################################################################################################################################################################################
"""
This function calculates soil metric potential from soil retention curve type and soil volumetric water potential. The supported methods are

$(METHODLIST)
"""
function soil_ψ_25 end


"""
    soil_ψ_25(vc::Union{BrooksCorey{FT}, VanGenuchten{FT}}, θ::FT) where {FT<:AbstractFloat}

Return the soil metric potential, given
- `vc` [`BrooksCorey`](@ref) or [`VanGenuchten`](@ref) type structure
- `θ` Soil volumetric water content (absolute value)

# Examples
```
soil_ψ_25(bc, FT(0.2))
soil_ψ_25(vg, FT(0.2))
```
"""
soil_ψ_25(vc::BrooksCorey{FT}, θ::FT) where {FT<:AbstractFloat} = (
    @unpack B, Ψ_SAT, Θ_RES, Θ_SAT = vc;

    # calculate effective θ
    _θ_e = min(1, max(0, (θ - Θ_RES) / (Θ_SAT - Θ_RES)));

    # if _θ_e >= 1, return 0
    if _θ_e >= 1
        return FT(0)
    end;

    return -Ψ_SAT / (_θ_e ^ B)
)


soil_ψ_25(vc::VanGenuchten{FT}, θ::FT) where {FT<:AbstractFloat} = (
    @unpack M, N, Α, Θ_RES, Θ_SAT = vc;

    # calculate effective θ
    _θ_e = min(1, max(0, (θ - Θ_RES) / (Θ_SAT - Θ_RES)));

    # if _θ_e >= 1, return 0
    if _θ_e >= 1
        return FT(0)
    end;

    return -1 * (_θ_e ^ (-1/M) - 1) ^ (1/N) / Α
)


end # module
