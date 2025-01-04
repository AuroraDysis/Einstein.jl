module FDMSuite

using ArgCheck: @argcheck
using LinearAlgebra

include("grid.jl")
include("weights.jl")
include("dissipation.jl")

end
