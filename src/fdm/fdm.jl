module FDMSuite

using ArgCheck: @argcheck
using LinearAlgebra

include("utils.jl")
include("grid.jl")
include("weights.jl")
include("stencil.jl")
include("diss.jl")

end
