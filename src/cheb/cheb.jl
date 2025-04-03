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

# Chebyshev series
include("series/evaluate.jl")
include("series/derivative.jl")
include("series/integrate.jl")
include("series/integration_matrix.jl")
include("series/chop.jl")

# Gauss-Chebyshev grid
include("chebyshev_gauss/grid.jl")

# Gauss-Chebyshev-Lobatto grid
include("chebyshev_lobatto/grid.jl")

include("dissipation.jl")

end
