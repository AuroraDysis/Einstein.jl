"""
    cheb_feval(f::VT, x::TR) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}

Evaluate a Chebyshev series at specified points.

# Mathematical Background
For a Chebyshev series:
```math
f(x) = \\sum_{k=0}^n c_k T_k(x)
```
where ``T_k(x)`` are Chebyshev polynomials of the 1st kind, this module provides
two evaluation methods:

1. Clenshaw algorithm (default) - Efficient for moderate degree polynomials
2. NDCT (TODO) - Fast for high degree polynomials with many evaluation points

# Current Implementation
Currently uses Clenshaw's algorithm exclusively. Future versions will implement NDCT
for improved performance on high-degree polynomials (n > 4000) or many evaluation points.

# Performance Notes
- Clenshaw's algorithm: O(n) operations per point
- NDCT (planned): O(n log n) operations for many points simultaneously

# Examples
```julia
# Evaluate using coefficients directly
coeffs = [1.0, 2.0, 3.0]  # f(x) = 1 + 2T₁(x) + 3T₂(x)
x = 0.5
y = cheb_clenshaw(coeffs, x)
```

# References
1. [doi:10.1137/1.9781611975949](@citet*)
2. [mason2002chebyshev](@citet*)
3. [chebfun/@chebtech/feval.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech/feval.m)
"""
function cheb_feval(f::VT, x::TR) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}
    return cheb_clenshaw(f, x)
end

export cheb_feval
