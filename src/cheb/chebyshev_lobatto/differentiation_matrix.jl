"""
    cheb_lobatto_differentiation_matrix([TR=Float64], n::Integer, k::Integer=1) where {TR<:AbstractFloat}
    cheb_lobatto_differentiation_matrix([TR=Float64], n::Integer, domain_width::TR, k::Integer=1) where {TR<:AbstractFloat}

Construct a Chebyshev differentiation that maps function values at `n` Chebyshev points of the 2nd kind 
to values of the `k`-th derivative of the interpolating polynomial at those points.

# Arguments
- `TR`: Element type (defaults to Float64)
- `n::Integer`: Number of Chebyshev points
- `k::Integer=1`: Order of the derivative (default: 1)

# References
- [chebfun/@chebcolloc2/chebcolloc2.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebcolloc2/chebcolloc2.m)
"""
function cheb_lobatto_differentiation_matrix(
    ::Type{TR}, n::Integer, k::Integer=1
) where {TR<:AbstractFloat}
    x = cheb_lobatto_points(TR, n)               # First kind points.
    w = cheb_lobatto_barycentric_weights(TR, n)           # Barycentric weights.
    t = cheb_lobatto_angles(TR, n)            # acos(x).
    D = barycentric_differentiation_matrix(x, w, k, t)       # Construct matrix.
    return D
end

function cheb_lobatto_differentiation_matrix(n::Integer, k::Integer=1)
    return cheb_lobatto_differentiation_matrix(Float64, n, k)
end

function cheb_lobatto_differentiation_matrix(
    ::Type{TR}, n::Integer, domain_width::TR, k::Integer=1
) where {TR<:AbstractFloat}
    D = cheb_lobatto_differentiation_matrix(TR, n, k)
    jacobian = (2 / domain_width)^k
    D .*= jacobian
    return D
end

function cheb_lobatto_differentiation_matrix(
    n::Integer, domain_width::Float64, k::Integer=1
)
    return cheb_lobatto_differentiation_matrix(Float64, n, domain_width, k)
end

export cheb_lobatto_differentiation_matrix
