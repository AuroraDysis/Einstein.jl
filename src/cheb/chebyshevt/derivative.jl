"""
    chebyshevt_derivative(coeffs::AbstractVector{T}) where {T<:AbstractFloat}
    chebyshevt_derivative!(coeffs::AbstractVector{T}, coeffs_der::AbstractVector{T}) where {T<:AbstractFloat}

Compute derivatives of Chebyshev coefficients.

# Arguments
- `coeffs`: Input vector of Chebyshev coefficients with length n
- `coeffs_der`: Pre-allocated output vector for derivative coefficients (length at least n - 1)
"""
function chebyshevt_derivative!(
    coeffs::AbstractVector{T}, coeffs_der::AbstractVector{T}
) where {T<:AbstractFloat}
    n = length(coeffs)
    n_der = length(coeffs_der)

    @argcheck n >= 1 "coeffs must have at least one element"
    @argcheck n_der >= n - 1 "coeffs_der must have at least n - 1 elements"

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
        half = one(T) / 2
        coeffs_der[1] *= half
    end

    return nothing
end

function chebyshevt_derivative(coeffs::AbstractVector{T}) where {T<:AbstractFloat}
    n = length(coeffs)
    coeffs_der = similar(coeffs, n - 1)
    chebyshevt_derivative!(coeffs, coeffs_der)
    return coeffs_der
end

export chebyshevt_derivative, chebyshevt_derivative!
