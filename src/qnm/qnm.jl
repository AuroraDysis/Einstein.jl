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

@enumx SchwPType ReggeWheeler Zerilli
export SchwPType

include("schw/reggepole.jl")
include("schw/cheb.jl")
include("schw/di.jl")

include("kerr/cf.jl")
include("kerr/cheb.jl")
include("kerr/di.jl")

end
