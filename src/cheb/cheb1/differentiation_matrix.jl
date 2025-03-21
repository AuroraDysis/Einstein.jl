"""
    cheb1_differentiation_matrix([TR=Float64], n::Integer, k::Integer=1) where {TR<:AbstractFloat}

Construct a Chebyshev differentiation that maps function values at `n` Chebyshev points of the 1st kind 
to values of the `k`-th derivative of the interpolating polynomial at those points.

# Arguments
- `TR`: Element type (defaults to Float64)
- `n::Integer`: Number of Chebyshev points
- `k::Integer=1`: Order of the derivative (default: 1)

# References
- [chebfun/@chebcolloc1/chebcolloc1.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebcolloc1/chebcolloc1.m)
"""
function cheb1_differentiation_matrix(
    ::Type{TR}, n::Integer, k::Integer=1
) where {TR<:AbstractFloat}
    x = cheb1_points(TR, n)               # First kind points.
    w = cheb1_barycentric_weights(TR, n)           # Barycentric weights.
    t = cheb1_angles(TR, n)            # acos(x).
    D = barycentric_differentiation_matrix(x, w, k, t)       # Construct matrix.
    return D
end

function cheb1_differentiation_matrix(n::Integer, k::Integer=1)
    return cheb1_differentiation_matrix(Float64, n, k)
end

function cheb1_differentiation_matrix(
    ::Type{TR}, n::Integer, lower_bound::TR, upper_bound::TR, k::Integer=1
) where {TR<:AbstractFloat}
    D = cheb1_differentiation_matrix(TR, n, k)
    scale = (2 / (upper_bound - lower_bound))^k
    D .*= scale
    return D
end

function cheb1_differentiation_matrix(
    n::Integer, lower_bound::Float64, upper_bound::Float64, k::Integer=1
)
    return cheb1_differentiation_matrix(Float64, n, lower_bound, upper_bound, k)
end

function _cheb_differentiation_matrix(
    ::ChebyshevFirstKindNode, ::Type{TR}, n::Integer, lower_bound::TR, upper_bound::TR, k::Integer
) where {TR<:AbstractFloat}
    return cheb1_differentiation_matrix(TR, n, lower_bound, upper_bound, k)
end

export cheb1_differentiation_matrix
