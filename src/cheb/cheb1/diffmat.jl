"""
    cheb1_diffmat([TR=Float64], n::Integer, k::Integer=1) where {TR<:AbstractFloat}

Construct a Chebyshev differentiation that maps function values at `n` Chebyshev points of the 1st kind 
to values of the `k`-th derivative of the interpolating polynomial at those points.

# Arguments
- `TR`: Element type (defaults to Float64)
- `n::Integer`: Number of Chebyshev points
- `k::Integer=1`: Order of the derivative (default: 1)

# References
- [chebfun/@chebcolloc1/chebcolloc1.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebcolloc1/chebcolloc1.m)
"""
function cheb1_diffmat(::Type{TR}, n::Integer, k::Integer=1) where {TR<:AbstractFloat}
    x = cheb1_pts(TR, n)               # First kind points.
    w = cheb1_barywts(TR, n)           # Barycentric weights.
    t = cheb1_angles(TR, n)            # acos(x).
    D = bary_diffmat(x, w, k, t)       # Construct matrix.
    return D
end

function cheb1_diffmat(n::Integer, k::Integer=1)
    return cheb1_diffmat(Float64, n, k)
end

export cheb1_diffmat
