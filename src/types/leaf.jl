"""
$(TYPEDEF)

Struct that contains leaf traits such as plant photosynthesis, hydraulics, and biophysics.

# Fields
$(TYPEDFIELDS)
"""
mutable struct Leaf{FT<:AbstractFloat}
    # variables that changes with time
    "Current leaf temperature"
    t::FT

    # variables that are used as a cache
    "Last leaf temperature. If different from t, then make temperature correction"
    _t_old::FT

    # constructor
    Leaf{FT}() where {FT<:AbstractFloat} = new{FT}(298.15, 0.0);
end
