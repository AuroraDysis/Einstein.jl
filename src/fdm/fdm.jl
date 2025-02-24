module FiniteDifferenceSuite

using ..Utils

using ArgCheck: @argcheck

using EnumX
using FillArrays
using FastBroadcast
using LinearAlgebra
using BandedMatrices
using StaticArrays

include("utils.jl")
include("grid.jl")
include("weights.jl")
include("dissipation.jl")
include("operator.jl")
include("integrate.jl")
include("interpolation.jl")

end
