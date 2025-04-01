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

include("grid.jl")
include("angles.jl")
include("interpolation.jl")
include("analysis.jl")
include("synthesis.jl")
include("analysis_matrix.jl")
include("synthesis_matrix.jl")
include("differentiation_matrix.jl")
include("integration_matrix.jl")
include("dissipation.jl")
include("quadrature_weights.jl")

include("utils/fft.jl")
include("utils/barycentric_differentiation_matrix.jl")
include("utils/coeffs_diff.jl")
include("utils/coeffs_eval.jl")
include("utils/coeffs_integrate.jl")
include("utils/coeffs_integration_matrix.jl")

include("chebtech2/points.jl")
include("chebtech2/angles.jl")
include("chebtech2/synthesis.jl")
include("chebtech2/analysis.jl")
include("chebtech2/barycentric_weights.jl")
include("chebtech2/quadrature_weights.jl")
include("chebtech2/differentiation_matrix.jl")
include("chebtech2/analysis_matrix.jl")
include("chebtech2/synthesis_matrix.jl")
include("chebtech2/integration_matrix.jl")

# Rectangular spectral collocation
include("rect/differentiation_matrix.jl")
include("rect/integration_matrix.jl")

# Ultraspherical spectral method
include("ultra/diffmat.jl")
include("ultra/convertmat.jl")
include("ultra/sphankel.jl")
include("ultra/multmat.jl")

end
