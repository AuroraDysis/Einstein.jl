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

    # Moments - Exact integrals of T_k (even)
    nm = div(n - 1, 2) + 1
    m = Vector{TF}(undef, nm)
    @inbounds begin
        m[1] = 2
        for i in 2:nm
            m[i] = 2 / (one(TF) - (2 * (i - 1))^2)
        end
    end

    # Mirror the vector for the use of ifft
    c = Vector{Complex{TF}}(undef, n)

    c[1:nm] .= m
    if isodd(n)
        start_idx = div(n + 1, 2)
        @inbounds for i in (nm + 1):n
            c[i] = -m[start_idx - i + nm + 1]
        end
    else
        start_idx = div(n, 2)
        c[nm + 1] = 0
        @inbounds for i in (nm + 2):n
            c[i] = -m[start_idx - i + nm + 2]
        end
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
    ::Type{TF}, n::Integer, lower_bound::TF, upper_bound::TF
) where {TF<:AbstractFloat}
    weights = cheb_gauss_quadrature_weights(TF, n)
    jacobian = (upper_bound - lower_bound) / 2
    weights .*= jacobian
    return weights
end

function cheb_gauss_quadrature_weights(
    n::Integer, lower_bound::Float64, upper_bound::Float64
)
    return cheb_gauss_quadrature_weights(Float64, n, lower_bound, upper_bound)
end

export cheb_gauss_quadrature_weights
