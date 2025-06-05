"""
    cheb_gauss_quadrature_weights([TF=Float64], n::Integer) where {TF<:AbstractFloat}

Compute quadrature weights for Chebyshev points of the 1st kind.

# Arguments
- `TF`: Type parameter for the weights (e.g., Float64)
- `n`: Number of points

# References
- [chebfun/@chebtech1/quadwts.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech1/quadwts.m)
"""
function cheb_gauss_quadrature_weights(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    if n == 0
        return TF[]
    elseif n == 1
        return TF[2]
    end

    c = Vector{Complex{TF}}(undef, n)

    # Moments - Exact integrals of T_k (even)
    nm = div(n - 1, 2) + 1
    c[1] = 2
    @inbounds for i in 2:nm
        c[i] = 2 / (one(TF) - (2 * (i - 1))^2)
    end

    # Mirror the vector for the use of ifft
    if isodd(n)
        start_idx = div(n + 1, 2)
        @.. c[(nm + 1):n] = -@view(c[start_idx:-1:2])
    else
        c[nm + 1] = 0
        start_idx = div(n, 2)
        @.. c[(nm + 2):n] = -@view(c[start_idx:-1:2])
    end

    # Apply weight (rotation) vector
    im_pi_over_n = im * convert(TF, π) / n
    @inbounds for k in 0:(n - 1)
        c[k + 1] *= exp(k * im_pi_over_n)
    end

    ifft!(c)
    weights = real.(c)

    return weights
end

function cheb_gauss_quadrature_weights(n::Integer)
    return cheb_gauss_quadrature_weights(Float64, n)
end

function cheb_gauss_quadrature_weights(
    ::Type{TF}, n::Integer, domain_width::TF
) where {TF<:AbstractFloat}
    weights = cheb_gauss_quadrature_weights(TF, n)
    jacobian = domain_width / 2
    weights .*= jacobian
    return weights
end

function cheb_gauss_quadrature_weights(
    n::Integer, domain_width::Float64
)
    return cheb_gauss_quadrature_weights(Float64, n, domain_width)
end

export cheb_gauss_quadrature_weights
