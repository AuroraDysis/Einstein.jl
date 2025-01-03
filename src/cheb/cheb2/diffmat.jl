"""
    cheb2_diffmat([TR=Float64], n::TI, k::TI=1) where {TR<:AbstractFloat,TI<:Integer}

Construct a Chebyshev differentiation that maps function values at `n` Chebyshev points of the 2nd kind 
to values of the `k`-th derivative of the interpolating polynomial at those points.

# Arguments
- `TR`: Element type (defaults to Float64)
- `n::Integer`: Number of Chebyshev points
- `k::Integer=1`: Order of the derivative (default: 1)

# References
- [chebfun/@chebcolloc2/chebcolloc2.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebcolloc2/chebcolloc2.m)
"""
function cheb2_diffmat(::Type{TR}, n::TI, k::TI=1) where {TR<:AbstractFloat,TI<:Integer}
    x = cheb2_pts(TR, n)               # First kind points.
    w = cheb2_barywts(TR, n)           # Barycentric weights.
    t = cheb2_angles(TR, n)            # acos(x).
    D = bary_diffmat(x, w, k, t)       # Construct matrix.
    return D
end

function cheb2_diffmat(n::TI, k::TI=1) where {TI<:Integer}
    return cheb2_diffmat(Float64, n, k)
end

export cheb2_diffmat
