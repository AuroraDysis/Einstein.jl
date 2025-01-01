module ChebyshevSuite

using InlineTest
using ArgCheck: @argcheck
using FastBroadcast: @..
using LinearAlgebra: dot

include("utils.jl")

include("cheb1/grid.jl")
include("cheb1/bary.jl")

include("cheb2/grid.jl")
include("cheb2/bary.jl")
include("cheb2/quad.jl")
include("cheb2/asmat.jl")
include("cheb2/intmat.jl")

end
