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

include("utils/fft.jl")
include("utils/diff.jl")
include("utils/bary_diffmat.jl")
include("utils/coeffs_eval.jl")
include("utils/coeffs_integrate.jl")
include("utils/coeffs_integration_matrix.jl")
include("utils/dissipation.jl")

include("cheb1/points.jl")
include("cheb1/angles.jl")
include("cheb1/synthesis.jl")
include("cheb1/analysis.jl")
include("cheb1/barycentric_weights.jl")
include("cheb1/quadrature_weights.jl")
include("cheb1/differentiation_matrix.jl")
include("cheb1/analysis_matrix.jl")
include("cheb1/synthesis_matrix.jl")
include("cheb1/integration_matrix.jl")

include("cheb2/points.jl")
include("cheb2/angles.jl")
include("cheb2/synthesis.jl")
include("cheb2/analysis.jl")
include("cheb2/barycentric_weights.jl")
include("cheb2/quadrature_weights.jl")
include("cheb2/differentiation_matrix.jl")
include("cheb2/analysis_matrix.jl")
include("cheb2/synthesis_matrix.jl")
include("cheb2/integration_matrix.jl")

# Rectangular spectral collocation
include("rect/rectint.jl")
include("rect/rectdiff.jl")

# Ultraspherical spectral method
include("ultra/diffmat.jl")
include("ultra/convertmat.jl")
include("ultra/sphankel.jl")
include("ultra/multmat.jl")

end
