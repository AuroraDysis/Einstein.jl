"""
The `GaussChebyshev` module provides a comprehensive set of tools for working with Chebyshev points of the 1st kind
and their associated operations. This module is particularly useful for:

1. Function Approximation:
   - Converting between function values at Chebyshev points and Chebyshev series coefficients
   - Computing barycentric weights for efficient interpolation
   - Generating Chebyshev points and their corresponding angles

2. Numerical Operations:
   - Computing differentiation matrices for Chebyshev collocation
   - Computing integration matrices for Chebyshev quadrature
   - Computing quadrature weights for numerical integration

The module is designed to work with the standard interval [-1,1] by default, but also supports
mapped intervals [a,b] through appropriate scaling transformations.

# References
- [chebfun/chebfun](https://github.com/chebfun/chebfun)
"""
module GaussChebyshev
include("angles.jl")
include("points.jl")

include("coeffs2vals.jl")
include("vals2coeffs.jl")

include("barycentric_weights.jl")
include("quadrature_weights.jl")

include("differentiation_matrix.jl")
include("integration_matrix.jl")
end

export GaussChebyshev
