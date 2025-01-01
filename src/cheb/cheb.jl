module ChebyshevSuite

using InlineTest

using FFTW: fft!, ifft!
using FastTransforms: fft!, ifft!
using ArgCheck: @argcheck
using FastBroadcast: @..

include("utils.jl")

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
