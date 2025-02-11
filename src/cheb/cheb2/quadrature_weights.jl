"""
    cheb2_quadrature_weights([TF=Float64], n::Integer) where {TF<:AbstractFloat}

Compute quadrature weights for Chebyshev points of the 2nd kind.

# Arguments
- `TF`: Type parameter for the weights (e.g., Float64)
- `n`: Number of points

# References
- [chebfun/@chebtech2/quadwts.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/quadwts.m)
"""
function cheb2_quadrature_weights(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
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

function cheb2_quadrature_weights(n::Integer)
    return cheb2_quadrature_weights(Float64, n)
end

function cheb2_quadrature_weights(
    ::Type{TF}, n::Integer, x_min::TF, x_max::TF
) where {TF<:AbstractFloat}
    w = cheb2_quadrature_weights(TF, n)
    scale = (x_max - x_min) / 2
    w .*= scale
    return w
end

function cheb2_quadrature_weights(n::Integer, x_min::Float64, x_max::Float64)
    return cheb2_quadrature_weights(Float64, n, x_min, x_max)
end

function _cheb_quadrature_weights(
    ::ChebyshevSecondKindNode, ::Type{TF}, n::Integer, x_min::TF, x_max::TF
) where {TF<:AbstractFloat}
    return cheb2_quadrature_weights(TF, n, x_min, x_max)
end

export cheb2_quadrature_weights
