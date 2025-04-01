"""
The `GaussChebyshev` module provides a comprehensive set of tools for working with Chebyshev points of the 1st kind
and coefficients of the corresponding 1st-kind Chebyshev series expansion.

The module is designed to work with the standard interval [-1,1] by default, but also supports
mapped intervals [a,b] through appropriate scaling transformations.

# References
- [chebfun/@chebtech1/chebtech1.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech1/chebtech1.m)
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
