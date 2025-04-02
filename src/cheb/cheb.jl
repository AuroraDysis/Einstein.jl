module ChebyshevSuite

using AbstractFFTs: Plan
using ArgCheck: @argcheck
using FastBroadcast: @..

using FFTW
using LinearAlgebra
using FillArrays
using SparseArrays
using BandedMatrices
using ToeplitzMatrices: Toeplitz, Hankel

using ..Utils

include("utils/fft.jl")
include("utils/barycentric_differentiation_matrix.jl")

# utils for coefficients of the corresponding first-kind Chebyshev series expansion.
include("chebyshevt/evaluate.jl")
include("chebyshevt/derivative.jl")
include("chebyshevt/integrate.jl")
include("chebyshevt/integration_matrix.jl")

# Gauss-Chebyshev grid
include("chebyshev_gauss/grid.jl")

# Gauss-Chebyshev-Lobatto grid
include("chebyshev_lobatto/grid.jl")

# Rectangular spectral collocation
include("rect/differentiation_matrix.jl")
include("rect/integration_matrix.jl")

# Ultraspherical spectral method
include("ultra/diffmat.jl")
include("ultra/convertmat.jl")
include("ultra/sphankel.jl")
include("ultra/multmat.jl")

include("dissipation.jl")

end
