module QNMSuite

using ..Utils
using ..ChebyshevSuite

using ArgCheck: @argcheck
using SciMLBase: AbstractNonlinearAlgorithm, AbstractODEAlgorithm

using FillArrays
using StaticArrays
using LinearAlgebra
using FastBroadcast
using NonlinearSolve

include("utils/continued_fraction.jl")
include("utils/polyeig.jl")

include("schw/qnm.jl")
include("kerr/qnm.jl")

end
