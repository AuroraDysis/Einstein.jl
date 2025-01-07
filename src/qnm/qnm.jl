module QNMSuite

using ..Utils
using ..ChebSuite

using ArgCheck: @argcheck
using SciMLBase: AbstractNonlinearAlgorithm, AbstractODEAlgorithm

using EnumX
using ApproxFun
using FillArrays
using Parameters
using StaticArrays
using LinearAlgebra
using FastBroadcast
using NonlinearSolve
using OrdinaryDiffEq

@enumx BCType Natural Dirichlet

export BCType

include("utils/contfrac.jl")
include("utils/polyeig.jl")

include("schw/schw.jl")

include("kerr/kerr.jl")

end
