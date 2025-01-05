module QNMSuite

using ..Utils
using ..ChebSuite

using ArgCheck: @argcheck
using NonlinearEigenproblems: PEP

using EnumX
using ApproxFun
using Parameters
using LinearAlgebra

include("utils/contfrac.jl")

include("schw/ultra.jl")

end
