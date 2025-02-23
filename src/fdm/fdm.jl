module FiniteDifferenceSuite

using ..Utils

using ArgCheck: @argcheck

using FillArrays
using FastBroadcast
using LinearAlgebra
using BandedMatrices

include("utils.jl")
include("grid.jl")
include("weights.jl")
include("dissipation.jl")
include("operator.jl")
include("differentiation_matrix.jl")
include("integrate.jl")
include("interpolation.jl")

end
