module QNMSuite

using ..Utils

using ArgCheck: @argcheck
using SciMLBase: AbstractNonlinearAlgorithm, AbstractODEAlgorithm

using FillArrays
using StaticArrays
using SparseArrays
using LinearAlgebra
using FastBroadcast
using NonlinearSolve

include("utils/continued_fraction.jl")
include("utils/companion.jl")

include("schw/qnm.jl")
include("kerr/qnm.jl")

end
