module ChebyshevSuite

using InlineTest
using ArgCheck: @argcheck
using FastBroadcast: @..
using LinearAlgebra: dot

include("utils.jl")

include("cheb1/angles.jl")
include("cheb1/chebpts.jl")
include("cheb1/bary.jl")
include("cheb1/quad.jl")

include("cheb2/angles.jl")
include("cheb2/chebpts.jl")
include("cheb2/bary.jl")
include("cheb2/quad.jl")
include("cheb2/asmat.jl")
include("cheb2/intmat.jl")

end
