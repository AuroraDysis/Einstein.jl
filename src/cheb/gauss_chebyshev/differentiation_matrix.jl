"""
    gauss_chebyshev_differentiation_matrix([TR=Float64], n::Integer, k::Integer=1) where {TR<:AbstractFloat}

Construct a Chebyshev differentiation that maps function values at `n` Chebyshev points of the 1st kind 
to values of the `k`-th derivative of the interpolating polynomial at those points.

# Arguments
- `TR`: Element type (defaults to Float64)
- `n::Integer`: Number of Chebyshev points
- `k::Integer=1`: Order of the derivative (default: 1)

# References
- [chebfun/@chebcolloc1/chebcolloc1.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebcolloc1/chebcolloc1.m)
"""
function gauss_chebyshev_differentiation_matrix(
    ::Type{TR}, n::Integer, k::Integer=1
) where {TR<:AbstractFloat}
    x = gauss_chebyshev_points(TR, n)
    w = gauss_chebyshev_barycentric_weights(TR, n)
    t = gauss_chebyshev_angles(TR, n)
    D = barycentric_differentiation_matrix(x, w, k, t)
    return D
end

function gauss_chebyshev_differentiation_matrix(n::Integer, k::Integer=1)
    return gauss_chebyshev_differentiation_matrix(Float64, n, k)
end

function gauss_chebyshev_differentiation_matrix(
    ::Type{TR}, n::Integer, lower_bound::TR, upper_bound::TR, k::Integer=1
) where {TR<:AbstractFloat}
    D = gauss_chebyshev_differentiation_matrix(TR, n, k)
    jacobian = (2 / (upper_bound - lower_bound))^k
    D .*= jacobian
    return D
end

function gauss_chebyshev_differentiation_matrix(
    n::Integer, lower_bound::Float64, upper_bound::Float64, k::Integer=1
)
    return gauss_chebyshev_differentiation_matrix(Float64, n, lower_bound, upper_bound, k)
end

export gauss_chebyshev_differentiation_matrix
