module ChebSuite

import AbstractFFTs: Plan
using ..Utils
using FFTW
using FastTransforms
using ArgCheck: @argcheck
using FastBroadcast: @..
using LinearAlgebra
using FillArrays
using SparseArrays
using BandedMatrices
using ToeplitzMatrices: Toeplitz, Hankel

include("utils/bary.jl")
include("utils/clenshaw.jl")
include("utils/feval.jl")
include("utils/fft.jl")
include("utils/diff.jl")
include("utils/cumsum.jl")
include("utils/bary_diffmat.jl")
include("utils/cumsummat.jl")
include("utils/dissipation.jl")

include("cheb1/grid.jl")
include("cheb1/angles.jl")
include("cheb1/pts.jl")
include("cheb1/coeffs2vals.jl")
include("cheb1/vals2coeffs.jl")
include("cheb1/barywts.jl")
include("cheb1/interp.jl")
include("cheb1/quad.jl")
include("cheb1/diffmat.jl")
include("cheb1/asmat.jl")
include("cheb1/cumsummat.jl")

include("cheb1/grid.jl")
include("cheb2/angles.jl")
include("cheb2/pts.jl")
include("cheb2/coeffs2vals.jl")
include("cheb2/vals2coeffs.jl")
include("cheb2/barywts.jl")
include("cheb2/interp.jl")
include("cheb2/quad.jl")
include("cheb2/diffmat.jl")
include("cheb2/asmat.jl")
include("cheb2/cumsummat.jl")

# Rectangular spectral collocation
include("rect/rectint.jl")
include("rect/rectdiff.jl")

# Ultraspherical spectral method
include("ultra/diffmat.jl")
include("ultra/convertmat.jl")
include("ultra/sphankel.jl")
include("ultra/multmat.jl")

end
