"""
    cheb_lobatto_quadrature_weights([TF=Float64], n::Integer) where {TF<:AbstractFloat}

Compute quadrature weights for Chebyshev points of the 2nd kind.

# Arguments
- `TF`: Type parameter for the weights (e.g., Float64)
- `n`: Number of points

# References
- [chebfun/@chebtech2/quadwts.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/quadwts.m)
"""
function cheb_lobatto_quadrature_weights(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    if n == 0
        return TF[]
    elseif n == 1
        return TF[2]
    end

    c = Array{Complex{TF}}(undef, n - 1)

    # Moments - Exact integrals of T_k (even)
    nm = div(n - 1, 2) + 1
    c[1] = 2
    @inbounds for i in 2:nm
        c[i] = 2 / (one(TF) - (2 * (i - 1))^2)
    end

    # Mirror for DCT via FFT
    half_idx = div(n, 2)
    c[(nm + 1):(n - 1)] .= @view(c[half_idx:-1:2])

    ifft!(c)

    weights = Array{TF}(undef, n)
    @.. weights[1:(n - 1)] = real(c)
    bc_weight = weights[1] / 2
    weights[1] = bc_weight
    weights[n] = bc_weight

    return weights
end

function cheb_lobatto_quadrature_weights(n::Integer)
    return cheb_lobatto_quadrature_weights(Float64, n)
end

function cheb_lobatto_quadrature_weights(
    ::Type{TF}, n::Integer, lower_bound::TF, upper_bound::TF
) where {TF<:AbstractFloat}
    weights = cheb_lobatto_quadrature_weights(TF, n)
    jacobian = (upper_bound - lower_bound) / 2
    weights .*= jacobian
    return weights
end

function cheb_lobatto_quadrature_weights(
    n::Integer, lower_bound::Float64, upper_bound::Float64
)
    return cheb_lobatto_quadrature_weights(Float64, n, lower_bound, upper_bound)
end

export cheb_lobatto_quadrature_weights
