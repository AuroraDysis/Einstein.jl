"""
    cheb_rectint([TR=Float64], n::Integer) where {TR<:AbstractFloat}
    cheb_rectint([TR=Float64], n::Integer, x_min::TR, x_max::TR) where {TR<:AbstractFloat}

Generate the Chebyshev integration matrix that operates directly on function values.

# Arguments
- `TR`: Type parameter for the matrix elements (e.g., Float64)
- `n`: Size of the matrix (n×n)
- `x_min`: (Optional) Lower bound of the integration interval
- `x_max`: (Optional) Upper bound of the integration interval

# References
- [chebfun/intmat.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/intmat.m)
"""
function cheb_rectint(::Type{TR}, m::Integer, n::Integer) where {TR<:AbstractFloat}
    # Build Lagrange basis
    K = Array{TR}(undef, n + 1, n)
    vals2coeffs_op = ChebyshevSecondKindAnalysis{TR}(n)
    @inbounds for i in 1:n
        K[1:(end - 1), i] = vals2coeffs_op(OneElement(one(TR), i, n))
    end

    # Integrate
    cumsum_op = ChebyshevCoefficientsCumulativeIntegration{TR}(n)
    @inbounds for i in 1:n
        K[:, i] = cumsum_op(@view(K[1:(end - 1), i]))
    end

    # Evaluate at grid
    xm = cheb2_points(m)
    intmat = Array{TR}(undef, m, n)
    @inbounds for j in 1:n, i in 1:n
        intmat[i, j] = cheb_feval(@view(K[:, j]), xm[i])
    end

    return intmat
end

function cheb_rectint(n::Integer)
    return cheb_rectint(Float64, n, n)
end

function cheb_rectint(
    ::Type{TR}, m::Integer, n::Integer, x_min::TR, x_max::TR
) where {TR<:AbstractFloat}
    intmat = cheb_rectint(TR, m, n)
    scale = (x_max - x_min) / 2
    intmat .*= scale
    return intmat
end

function cheb_rectint(m::Integer, n::Integer, x_min::Float64, x_max::Float64)
    return cheb_rectint(Float64, m, n, x_min, x_max)
end

function cheb_rectint(n::Integer, x_min::Float64, x_max::Float64)
    return cheb_rectint(Float64, n, n, x_min, x_max)
end

export cheb_rectint
