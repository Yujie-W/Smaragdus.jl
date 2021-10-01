"""
$(TYPEDEF)

Hierachy of AbstractSoilVC:
- [`BrooksCorey`](@ref)
- [`VanGenuchten`](@ref)
"""
abstract type AbstractSoilVC{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2021-Sep-30: define this structure with no default constructor
#     2021-Sep-30: define a constructor to fil the value from VanGenuchten parameters
#
#######################################################################################################################################################################################################
"""
$(TYPEDEF)

Brooks Corey soil parameters

# Fields
$(TYPEDFIELDS)

# Examples
```julia
bc = BrooksCorey{FT}("Test", FT(5), FT(2), FT(0.5), FT(0.1));
bc = BrooksCorey{FT}(VanGenuchten{FT}("Loam"));
```
"""
struct BrooksCorey{FT<:AbstractFloat} <:AbstractSoilVC{FT}
    "Soil type"
    TYPE::String
    "Soil b"
    B::FT
    "Potential at saturation `[MPa]`"
    Ψ_SAT::FT
    "Saturated soil volumetric water content"
    Θ_SAT::FT
    "Residual soil volumetric water content"
    Θ_RES::FT
end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2021-Sep-30: define this structure with two default constructors from an incomplete parameter set
#
#######################################################################################################################################################################################################
"""
$(TYPEDEF)

van Genuchten soil parameters

# Fields
$(TYPEDFIELDS)

# Examples
```julia
vg = VanGenuchten{FT}("Loam");
vg = VanGenuchten{FT}("Silt");
vg = VanGenuchten{FT}("Test", FT(100), FT(2), FT(0.5), FT(0.1));
```
"""
struct VanGenuchten{FT<:AbstractFloat} <:AbstractSoilVC{FT}
    "Soil type"
    TYPE::String
    "Soil α is related to the inverse of the air entry suction, α > 0"
    Α::FT
    "Soil n is Measure of the pore-size distribution"
    N::FT
    "Soil m = 1 - 1/n"
    M::FT
    "Saturated soil volumetric water content"
    Θ_SAT::FT
    "Residual soil volumetric water content"
    Θ_RES::FT

    # constructors
    VanGenuchten{FT}(name::String, α::FT, n::FT, θ_sat::FT, θ_res::FT) where {FT<:AbstractFloat} = new{FT}(name, α, n, 1-1/n, θ_sat, θ_res)

    VanGenuchten{FT}(name::String) where {FT<:AbstractFloat} = (
        # Parameters from Silt soil
        _p = [ 163.2656, 1.37, 0.46, 0.034];

        # switch name
        if name=="Sand"
            _p = [1479.5945, 2.68, 0.43, 0.045];
        elseif name=="Loamy Sand"
            _p = [1265.3084, 2.28, 0.41, 0.057];
        elseif name=="Sandy Loam"
            _p = [ 765.3075, 1.89, 0.41, 0.065];
        elseif name=="Loam"
            _p = [ 367.3476, 1.56, 0.43, 0.078];
        elseif name=="Sandy Clay Loam"
            _p = [ 602.0419, 1.48, 0.39, 0.100];
        elseif name=="Silt Loam"
            _p = [ 204.0820, 1.41, 0.45, 0.067];
        elseif name=="Silt"
            _p = [ 163.2656, 1.37, 0.46, 0.034];
        elseif name=="Clay Loam"
            _p = [ 193.8779, 1.31, 0.41, 0.095];
        elseif name=="Silty Clay Loam"
            _p = [ 102.0410, 1.23, 0.43, 0.089];
        elseif name== "Sandy Clay"
            _p = [ 275.5107, 1.23, 0.38, 0.100];
        elseif name=="Silty Clay"
            _p = [  51.0205, 1.09, 0.36, 0.070];
        elseif name=="Clay"
            _p = [  81.6328, 1.09, 0.38, 0.068];
        else
            @warn "Soil type $(name) not recognized, use Silt instead.";
            name = "Silt";
        end;

        # return a new struct
        return new{FT}(name, _p[1], _p[2], 1-1/_p[2], _p[3], _p[4])
    );
end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2021-Sep-30: define this structure with no constructors
#
#######################################################################################################################################################################################################
"""
$(TYPEDEF)

Struct that contains environmental conditions, such as soil moisture and atmospheric vapor pressure. Note that this structure is designed to be containers to interact with other CliMA modules and to
    prescribe values.

# Fields
$(TYPEDFIELDS)

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
    "Soil moisture retention curve"
    VC::Union{Vector{BrooksCorey{FT}}, Vector{VanGenuchten{FT}}}

    # variables that change with time

    # caches to speed up calculations

    # constructors
end
