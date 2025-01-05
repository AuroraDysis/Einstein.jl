module QNMSuite

using ..Utils
using ..ChebSuite

using ArgCheck: @argcheck
using NonlinearEigenproblems: PEP
using SciMLBase: AbstractNonlinearAlgorithm

using EnumX
using ApproxFun
using Parameters
using StaticArrays
using LinearAlgebra
using FastBroadcast
using NonlinearSolve

@enumx BCType Natural Dirichlet

export BCType

include("utils/contfrac.jl")

include("schw/schw.jl")

include("kerr/kerr.jl")

end
