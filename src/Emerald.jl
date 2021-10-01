#######################################################################################################################################################################################################
#
# Style Guide
# 1. Types and modules start with upper case letters for each word, such as LeafBiophysics
# 2. Functions are lower case words connected using _, such as leaf_biophysics
# 3. Constants are defined using all upper cases, such as LEAF_BIOPHYSICS
# 4. Variables are defined using all lower cases, such as leaf_bio_para
# 5. Temporary variables are defined to start with _, such as _leaf
# 6. Maximum length of each line is 200 letters (including space)
# 7. There should be 2 lines of  empty lines between different components, such as two functions and methods
# 8. Bug fixes or new changes should be documented in the comments above the struct, function, or method, such as this Style Guide above Emerald.jl
# 9. Function parameter list that spans multiple lines need to be spaced with 12 spaces (3 tabs)
#
#######################################################################################################################################################################################################
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
