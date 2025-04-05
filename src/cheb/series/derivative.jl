"""
    cheb_series_derivative!(coeffs::AbstractVector{TFC}) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    cheb_series_derivative!(coeffs_der::AbstractVector{TFC}, coeffs::AbstractVector{TFC}) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}

Compute derivatives of coefficients of Chebyshev series.

# Arguments
- `coeffs`: Input vector of Chebyshev coefficients with length n

# References
- [chebfun/@chebtech/diff.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech/diff.m)
"""
function cheb_series_derivative!(
    coeffs_der::AbstractVector{TFC}, coeffs::AbstractVector{TFC}
) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    @boundscheck begin
        @argcheck length(coeffs) >= 1 "coeffs must have at least one element"
        @argcheck length(coeffs_der) >= length(coeffs) - 1 "coeffs_der must have at least length(coeffs) - 1 elements"
    end

    n = length(coeffs)
    n_der = length(coeffs_der)

    if n == 1
        coeffs_der .= 0
        return nothing
    end

    # Compute 2k * c[k+1] for k = 1:n-1
    @inbounds for i in 1:(n - 1)
        coeffs_der[i] = 2 * i * coeffs[i + 1]
    end
    coeffs_der[n:end] .= 0

    # Compute cumulative sums for odd and even indices separately
    # This avoids the need to create temporary arrays
    @inbounds begin
        # Odd indices (n-1:-2:1)
        for i in (n - 1):-2:3
            coeffs_der[i - 2] += coeffs_der[i]
        end

        # Even indices (n-2:-2:1)
        for i in (n - 2):-2:3
            coeffs_der[i - 2] += coeffs_der[i]
        end

        # Scale first coefficient
        coeffs_der[1] /= 2
    end

    return nothing
end

function cheb_series_derivative!(
    coeffs::AbstractVector{TFC}
) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    cheb_series_derivative!(coeffs, coeffs)
    return nothing
end

export cheb_series_derivative!
