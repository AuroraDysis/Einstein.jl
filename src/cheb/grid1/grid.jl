"""
Module has two functions

1. approximating smooth functions on the interval [-1,1] using
function values at 1st-kind Chebyshev points and coefficients of the
corresponding 1st-kind Chebyshev series expansion.

2. Collocation discretization on 1st kind Chebyshev points.
"""
module GaussChebyshevGridSuite
include("angles.jl")
include("points.jl")

include("coeffs2vals.jl")
include("vals2coeffs.jl")

include("barycentric_weights.jl")
include("quadrature_weights.jl")

include("differentiation_matrix.jl")
include("integration_matrix.jl")
end

export GaussChebyshevGridSuite
