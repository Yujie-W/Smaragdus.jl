module Emerald

using Artifacts: @artifact_str
using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using MAT: matread
using Statistics: mean
using UnPack: @unpack


# export public types
export BrooksCorey, Leaf, VanGenuchten, WaveLengthSet


# include utility files
include("utils/radiation.jl" )
include("utils/statistics.jl")

# include type files
include("types/environment.jl")
include("types/wavelength.jl" )
include("types/leaf.jl"       )

# include sub-modules
include("Folium.jl")
include("Soli.jl"  )


end # module
