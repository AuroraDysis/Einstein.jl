"""
    cheb1_quadrature_weights([TR=Float64], n::Integer) where {TR<:AbstractFloat}

Compute quadrature weights for Chebyshev points of the 1st kind.

# Arguments
- `TR`: Type parameter for the weights (e.g., Float64)
- `n`: Number of points

# References
- [chebfun/@chebtech1/quadwts.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech1/quadwts.m)
"""
function cheb1_quadrature_weights(::Type{TR}, n::Integer) where {TR<:AbstractFloat}
    # Handle the special cases:
    if n == 0
        return TR[]
    elseif n == 1
        return TR[2]
    end

    # Preallocate the array m for the moments
    evens = 2:2:(n - 1)
    nm = length(evens) + 1
    m = Vector{TR}(undef, nm)
    @inbounds begin
        m[1] = 2 / one(TR)  # Corresponds to k=0
        for (i, k) in enumerate(evens)
            # m[i+1] = 1 - k^2
            m[i + 1] = 2 / (one(TR) - k^2)
        end
    end

    # Preallocate the coefficient array for FFT
    c = Vector{Complex{TR}}(undef, n)
    c[1:nm] .= m

    # Fill remaining coefficients based on parity of n
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

    # Multiply by rotation factors exp(1im*(0:n-1)*π/n)
    im_pi_over_n = im * convert(TR, π) / n
    @inbounds for k in 0:(n - 1)
        c[k + 1] *= exp(k * im_pi_over_n)
    end

    # Compute inverse FFT in-place and take the real part for the weights
    ifft!(c)
    w = real.(c)

    return w
end

function cheb1_quadrature_weights(n::Integer)
    return cheb1_quadrature_weights(Float64, n)
end

function cheb1_quadrature_weights(
    ::Type{TR}, n::Integer, x_min::TR, x_max::TR
) where {TR<:AbstractFloat}
    w = cheb1_quadrature_weights(TR, n)
    w .*= (x_max - x_min) / 2
    return w
end

function cheb1_quadrature_weights(n::Integer, x_min::Float64, x_max::Float64)
    return cheb1_quadrature_weights(Float64, n, x_min, x_max)
end

export cheb1_quadrature_weights
