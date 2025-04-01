"""
    chebgrid1_quadrature_weights([TF=Float64], n::Integer) where {TF<:AbstractFloat}

Compute quadrature weights for Chebyshev points of the 1st kind.

# Arguments
- `TF`: Type parameter for the weights (e.g., Float64)
- `n`: Number of points

# References
- [chebfun/@chebtech1/quadwts.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech1/quadwts.m)
"""
function chebgrid1_quadrature_weights(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    # Handle the special cases:
    if n == 0
        return TF[]
    elseif n == 1
        return TF[2]
    end

    # Preallocate the array m for the moments
    evens = 2:2:(n - 1)
    nm = length(evens) + 1
    m = Vector{TF}(undef, nm)
    @inbounds begin
        m[1] = 2 / one(TF)  # Corresponds to k=0
        for (i, k) in enumerate(evens)
            # m[i+1] = 1 - k^2
            m[i + 1] = 2 / (one(TF) - k^2)
        end
    end

    # Preallocate the coefficient array for FFT
    c = Vector{Complex{TF}}(undef, n)
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
    im_pi_over_n = im * convert(TF, π) / n
    @inbounds for k in 0:(n - 1)
        c[k + 1] *= exp(k * im_pi_over_n)
    end

    # Compute inverse FFT in-place and take the real part for the weights
    ifft!(c)
    w = real.(c)

    return w
end

function chebgrid1_quadrature_weights(n::Integer)
    return chebgrid1_quadrature_weights(Float64, n)
end

function chebgrid1_quadrature_weights(
    ::Type{TF}, n::Integer, lower_bound::TF, upper_bound::TF
) where {TF<:AbstractFloat}
    w = chebgrid1_quadrature_weights(TF, n)
    w .*= (upper_bound - lower_bound) / 2
    return w
end

function chebgrid1_quadrature_weights(n::Integer, lower_bound::Float64, upper_bound::Float64)
    return chebgrid1_quadrature_weights(Float64, n, lower_bound, upper_bound)
end
