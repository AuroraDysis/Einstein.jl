module ChebyshevSuite

using InlineTest

import AbstractFFTs: Plan
using FFTW
using FastTransforms
using ArgCheck: @argcheck
using FastBroadcast: @..
using LinearAlgebra

include("utils/bary.jl")
include("utils/clenshaw.jl")
include("utils/ndct.jl")
include("utils/fft.jl")
include("utils/cumsum.jl")

include("cheb1/angles.jl")
include("cheb1/pts.jl")
include("cheb1/coeffs2vals.jl")
include("cheb1/vals2coeffs.jl")
include("cheb1/barywts.jl")
include("cheb1/quad.jl")

include("cheb2/angles.jl")
include("cheb2/pts.jl")
include("cheb2/coeffs2vals.jl")
include("cheb2/vals2coeffs.jl")
include("cheb2/barywts.jl")
include("cheb2/quad.jl")
include("cheb2/asmat.jl")
include("cheb2/intmat.jl")

end
