include("angles.jl")
include("points.jl")

include("coeffs2vals.jl")
include("vals2coeffs.jl")

include("barycentric_weights.jl")
include("quadrature_weights.jl")

include("differentiation_matrix.jl")
include("integration_matrix.jl")

"""
The `GaussChebyshevGrid` module provides a comprehensive set of tools for working with Chebyshev points of the 1st kind
and coefficients of the corresponding 1st-kind Chebyshev series expansion.

The module is designed to work with the standard interval [-1,1] by default, but also supports
mapped intervals [a,b] through appropriate scaling transformations.

# References
- [chebfun/@chebtech1/chebtech1.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech1/chebtech1.m)
"""
module GaussChebyshevGrid
using ..ChebyshevSuite

const angles = gauss_chebyshev_angles
const points = gauss_chebyshev_points

const coeffs2vals = gauss_chebyshev_coeffs2vals
const vals2coeffs = gauss_chebyshev_vals2coeffs

const coeffs2vals_matrix = gauss_chebyshev_coeffs2vals_matrix
const vals2coeffs_matrix = gauss_chebyshev_vals2coeffs_matrix

const barycentric_weights = gauss_chebyshev_barycentric_weights
const quadrature_weights = gauss_chebyshev_quadrature_weights

const differentiation_matrix = gauss_chebyshev_differentiation_matrix
const integration_matrix = gauss_chebyshev_integration_matrix

export angles,
    points,
    coeffs2vals,
    vals2coeffs,
    coeffs2vals_matrix,
    vals2coeffs_matrix,
    barycentric_weights,
    quadrature_weights,
    differentiation_matrix,
    integration_matrix

end

export GaussChebyshevGrid
