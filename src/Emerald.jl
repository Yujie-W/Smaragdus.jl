module Emerald

using Artifacts: @artifact_str
using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using MAT: matread
using Statistics: mean
using UnPack: @unpack


# export public types
export Leaf, WaveLengthSet


# include files
include("utils/statistics.jl")

include("types/radiation.jl")
include("types/leaf.jl")


end # module
