"""
    cheb_feval(f::AbstractVector{TR}, x::TR) where {TR<:AbstractFloat}

Evaluate Chebyshev coefficients at a point.

# Performance Notes
- Clenshaw's algorithm: O(n) operations per point
- (TODO) NDCT: O(n log n) operations for many points simultaneously

# References
- [chebfun/@chebtech/feval.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech/feval.m)
"""
function cheb_feval(f::AbstractVector{TR}, x::TR) where {TR<:AbstractFloat}
    return cheb_clenshaw(f, x)
end

export cheb_feval
