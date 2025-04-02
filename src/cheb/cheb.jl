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

# utils for coefficients of the corresponding first-kind Chebyshev series expansion.
include("chebyshevt/evaluate.jl")
include("chebyshevt/derivative.jl")
include("chebyshevt/integrate.jl")
include("chebyshevt/integration_matrix.jl")

# Gauss-Chebyshev grid
include("chebyshev_gauss/grid.jl")

# Gauss-Chebyshev-Lobatto grid
include("chebyshev_lobatto/grid.jl")

include("dissipation.jl")

end
