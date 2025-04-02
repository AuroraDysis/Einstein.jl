module Utils

using ArgCheck: @argcheck
using FastBroadcast: @..

using Xsum
using FFTW
using LinearAlgebra

include("type.jl")
include("grid.jl")
include("errors.jl")
include("sum.jl")
include("dot.jl")
include("sws.jl")
include("fft.jl")

include("barycentric_weights.jl")
include("barycentric_interpolation.jl")
include("barycentric_differentiation_matrix.jl")

end
