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
include("fdmop.jl")
include("diffmat.jl")
include("dissmat.jl")
include("integrate.jl")
include("interpolation.jl")

end
