"""
    cheb_gauss_differentiation_matrix([TR=Float64], n::Integer, k::Integer=1) where {TR<:AbstractFloat}

Construct a Chebyshev differentiation that maps function values at `n` Chebyshev points of the 1st kind 
to values of the `k`-th derivative of the interpolating polynomial at those points.

# Arguments
- `TR`: Element type (defaults to Float64)
- `n::Integer`: Number of Chebyshev points
- `k::Integer=1`: Order of the derivative (default: 1)

# References
- [chebfun/@chebcolloc1/chebcolloc1.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebcolloc1/chebcolloc1.m)
"""
function cheb_gauss_differentiation_matrix(
    ::Type{TR}, n::Integer, k::Integer=1
) where {TR<:AbstractFloat}
    x = cheb_gauss_points(TR, n)
    w = cheb_gauss_barycentric_weights(TR, n)
    t = cheb_gauss_angles(TR, n)
    D = barycentric_differentiation_matrix(x, w, k, t)
    return D
end

function cheb_gauss_differentiation_matrix(n::Integer, k::Integer=1)
    return cheb_gauss_differentiation_matrix(Float64, n, k)
end

function cheb_gauss_differentiation_matrix(
    ::Type{TR}, n::Integer, lower_bound::TR, upper_bound::TR, k::Integer=1
) where {TR<:AbstractFloat}
    D = cheb_gauss_differentiation_matrix(TR, n, k)
    jacobian = (2 / (upper_bound - lower_bound))^k
    D .*= jacobian
    return D
end

function cheb_gauss_differentiation_matrix(
    n::Integer, lower_bound::Float64, upper_bound::Float64, k::Integer=1
)
    return cheb_gauss_differentiation_matrix(Float64, n, lower_bound, upper_bound, k)
end

export cheb_gauss_differentiation_matrix
