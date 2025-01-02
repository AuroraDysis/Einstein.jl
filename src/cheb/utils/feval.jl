"""
    cheb_feval(f::VT, x::TR) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}

Evaluate Chebyshev coefficients at a point.

# Performance Notes
- Clenshaw's algorithm: O(n) operations per point
- (TODO) NDCT: O(n log n) operations for many points simultaneously

# References
- [chebfun/@chebtech/feval.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech/feval.m)
"""
function cheb_feval(f::VT, x::TR) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}
    return cheb_clenshaw(f, x)
end

export cheb_feval
