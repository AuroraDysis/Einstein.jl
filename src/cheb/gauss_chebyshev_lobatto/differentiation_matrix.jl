"""
    gauss_chebyshev_lobatto_differentiation_matrix([TR=Float64], n::Integer, k::Integer=1) where {TR<:AbstractFloat}

Construct a Chebyshev differentiation that maps function values at `n` Chebyshev points of the 2nd kind 
to values of the `k`-th derivative of the interpolating polynomial at those points.

# Arguments
- `TR`: Element type (defaults to Float64)
- `n::Integer`: Number of Chebyshev points
- `k::Integer=1`: Order of the derivative (default: 1)

# References
- [chebfun/@chebcolloc2/chebcolloc2.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebcolloc2/chebcolloc2.m)
"""
function gauss_chebyshev_lobatto_differentiation_matrix(
    ::Type{TR}, n::Integer, k::Integer=1
) where {TR<:AbstractFloat}
    x = gauss_chebyshev_lobatto_points(TR, n)               # First kind points.
    w = gauss_chebyshev_lobatto_barycentric_weights(TR, n)           # Barycentric weights.
    t = gauss_chebyshev_lobatto_angles(TR, n)            # acos(x).
    D = barycentric_differentiation_matrix(x, w, k, t)       # Construct matrix.
    return D
end

function gauss_chebyshev_lobatto_differentiation_matrix(n::Integer, k::Integer=1)
    return gauss_chebyshev_lobatto_differentiation_matrix(Float64, n, k)
end

function gauss_chebyshev_lobatto_differentiation_matrix(
    ::Type{TR}, n::Integer, lower_bound::TR, upper_bound::TR, k::Integer=1
) where {TR<:AbstractFloat}
    D = gauss_chebyshev_lobatto_differentiation_matrix(TR, n, k)
    jacobian = (2 / (upper_bound - lower_bound))^k
    D .*= jacobian
    return D
end

function gauss_chebyshev_lobatto_differentiation_matrix(
    n::Integer, lower_bound::Float64, upper_bound::Float64, k::Integer=1
)
    return gauss_chebyshev_lobatto_differentiation_matrix(Float64, n, lower_bound, upper_bound, k)
end

export gauss_chebyshev_lobatto_differentiation_matrix
