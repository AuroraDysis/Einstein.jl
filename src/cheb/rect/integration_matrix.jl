"""
    cheb_rect_integration_matrix([TR=Float64], m::Integer, n::Integer) -> Matrix{TR}
    cheb_rect_integration_matrix([TR=Float64], m::Integer, n::Integer, lower_bound::TR, upper_bound::TR) -> Matrix{TR}

returns the m x n first-order rectangular integration matrix which maps from an n-point Chebyshev grid of the second kind to an m-point Chebyshev grid of the same kind.

# Arguments
- `TR`: Type parameter for the matrix elements (e.g., Float64)
- `m`: Size of the matrix (number of rows)
- `n`: Size of the matrix (number of columns)
- `lower_bound`: (Optional) Lower bound of the integration interval
- `upper_bound`: (Optional) Upper bound of the integration interval

# References
- [chebfun/intmat.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/intmat.m)
"""
function cheb_rect_integration_matrix(
    ::Type{TR}, m::Integer, n::Integer
) where {TR<:AbstractFloat}
    # Build Lagrange basis
    K = Array{TR}(undef, n + 1, n)
    vals2coeffs_op = GaussChebyshevLobatto.vals2coeffs(TR, n)
    @inbounds for i in 1:n
        K[1:(end - 1), i] = vals2coeffs_op(OneElement(one(TR), i, n))
    end

    # Integrate
    cumsum_op = chebyshevt_integrate(TR, n)
    @inbounds for i in 1:n
        K[:, i] = cumsum_op(@view(K[1:(end - 1), i]))
    end

    # Evaluate at grid
    xm = GaussChebyshevLobatto.points(TR, m)
    intmat = Array{TR}(undef, m, n)
    @inbounds for j in 1:n, i in 1:n
        intmat[i, j] = chebyshevt_evaluate(@view(K[:, j]), xm[i])
    end

    return intmat
end

function cheb_rect_integration_matrix(n::Integer)
    return cheb_rect_integration_matrix(Float64, n, n)
end

function cheb_rect_integration_matrix(
    ::Type{TR}, m::Integer, n::Integer, lower_bound::TR, upper_bound::TR
) where {TR<:AbstractFloat}
    intmat = cheb_rect_integration_matrix(TR, m, n)
    scale = (upper_bound - lower_bound) / 2
    intmat .*= scale
    return intmat
end

function cheb_rect_integration_matrix(
    m::Integer, n::Integer, lower_bound::Float64, upper_bound::Float64
)
    return cheb_rect_integration_matrix(Float64, m, n, lower_bound, upper_bound)
end

function cheb_rect_integration_matrix(n::Integer, lower_bound::Float64, upper_bound::Float64)
    return cheb_rect_integration_matrix(Float64, n, n, lower_bound, upper_bound)
end

export cheb_rect_integration_matrix
