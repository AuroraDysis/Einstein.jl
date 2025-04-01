module ChebyshevSuite

import AbstractFFTs: Plan
using ..Utils
using FFTW
using EnumX
using FastTransforms
using ArgCheck: @argcheck
using FastBroadcast: @..
using LinearAlgebra
using FillArrays
using SparseArrays
using BandedMatrices
using ToeplitzMatrices: Toeplitz, Hankel

include("interpolation.jl")
include("dissipation.jl")

include("utils/fft.jl")
include("utils/barycentric_differentiation_matrix.jl")
include("utils/coeffs_diff.jl")
include("utils/coeffs_eval.jl")
include("utils/coeffs_integrate.jl")
include("utils/coeffs_integration_matrix.jl")

include("grid/gauss_chebyshev/grid.jl")
include("grid/gauss_chebyshev_lobatto/grid.jl")

# Rectangular spectral collocation
include("rect/differentiation_matrix.jl")
include("rect/integration_matrix.jl")

# Ultraspherical spectral method
include("ultra/diffmat.jl")
include("ultra/convertmat.jl")
include("ultra/sphankel.jl")
include("ultra/multmat.jl")

end
