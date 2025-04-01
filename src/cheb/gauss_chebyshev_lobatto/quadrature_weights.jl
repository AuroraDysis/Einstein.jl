"""
    gauss_chebyshev_lobatto_quadrature_weights([TF=Float64], n::Integer) where {TF<:AbstractFloat}

Compute quadrature weights for Chebyshev points of the 2nd kind.

# Arguments
- `TF`: Type parameter for the weights (e.g., Float64)
- `n`: Number of points

# References
- [chebfun/@chebtech2/quadwts.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/quadwts.m)
"""
function gauss_chebyshev_lobatto_quadrature_weights(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    if n == 0
        return TF[]
    elseif n == 1
        return TF[2]
    end

    nm1 = n - 1

    # Fill exact integrals of T_k (even k)
    c = Array{Complex{TF}}(undef, nm1)
    @inbounds begin
        c[1] = TF(2)  # k = 0 case
        for k in 2:2:nm1
            c[k ÷ 2 + 1] = 2 / (one(TF) - k^2)
        end

        # Mirror for DCT via FFT (in-place)
        half_idx = floor(Int, n / 2)
        for i in 2:half_idx
            c[n - i + 1] = c[i]
        end

        # Compute weights via inverse FFT (in-place)
        ifft!(c)
    end

    # Adjust boundary weights (in-place)
    w = Vector{TF}(undef, n)
    @inbounds begin
        w[1] = real(c[1]) / 2
        for i in 2:nm1
            w[i] = real(c[i])
        end
        w[n] = real(c[1]) / 2
    end

    return w
end

function gauss_chebyshev_lobatto_quadrature_weights(n::Integer)
    return gauss_chebyshev_lobatto_quadrature_weights(Float64, n)
end

function gauss_chebyshev_lobatto_quadrature_weights(
    ::Type{TF}, n::Integer, lower_bound::TF, upper_bound::TF
) where {TF<:AbstractFloat}
    weights = gauss_chebyshev_lobatto_quadrature_weights(TF, n)
    jacobian = (upper_bound - lower_bound) / 2
    weights .*= jacobian
    return weights
end

function gauss_chebyshev_lobatto_quadrature_weights(n::Integer, lower_bound::Float64, upper_bound::Float64)
    return gauss_chebyshev_lobatto_quadrature_weights(Float64, n, lower_bound, upper_bound)
end

export gauss_chebyshev_lobatto_quadrature_weights
