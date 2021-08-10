module Emerald

using Artifacts: @artifact_str
using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using MAT: matread
using Statistics: mean
using UnPack: @unpack


# export public types
export Leaf, WaveLengthSet


# include utility files
include("utils/radiation.jl")
include("utils/statistics.jl")

# include type files
include("types/wavelength.jl")
include("types/leaf.jl")

# include sub-modules
include("LeafLevel/LeafLevel.jl")


end # module
