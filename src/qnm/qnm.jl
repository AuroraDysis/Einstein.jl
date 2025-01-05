module QNMSuite

using ..Utils
using ..ChebSuite

using ArgCheck: @argcheck
using NonlinearEigenproblems: PEP

using EnumX
using ApproxFun
using Parameters
using StaticArrays
using LinearAlgebra
using FastBroadcast

@enumx BCType Natural Dirichlet

export BCType

include("utils/contfrac.jl")

include("schw/expansion.jl")
include("schw/cheb.jl")

include("kerr/radial.jl")
include("kerr/cheb.jl")

end
