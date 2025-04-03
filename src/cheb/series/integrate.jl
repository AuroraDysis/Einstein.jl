@doc raw"""
    cheb_series_integrate(coeffs::AbstractVector{TF}) where {TF<:AbstractFloat}
    cheb_series_integrate(::Type{TF}, n::Integer)(coeffs::AbstractVector{TF}) where {TF<:AbstractFloat}

Compute the indefinite integral of a functio $f$ given its Chebyshev series,
with the constant of integration chosen such that $f(-1) = 0$.

# Arguments
- `coeffs`: Vector of Chebyshev series

# References
- [chebfun/@chebtech/cumsum.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech/cumsum.m)
"""
struct ChebyshevSeriesIntegrationPlan{TF<:AbstractFloat,TI<:Integer}
    n::TI              # Number of coefficients
    tmp::Vector{TF}    # Temporary storage for padded coefficients
    result::Vector{TF} # Result storage
    weights::Vector{TI}      # Pre-computed alternating signs

    function ChebyshevSeriesIntegrationPlan{TF}(n::TI) where {TF<:AbstractFloat,TI<:Integer}
        # Pre-allocate workspace
        tmp = Vector{TF}(undef, n + 2)
        result = Vector{TF}(undef, n + 1)

        # Pre-compute alternating signs [1, -1, 1, -1, ...]
        weights = ones(TI, n)
        weights[2:2:end] .= -one(TI)

        return new{TF,TI}(n, tmp, result, weights)
    end
end

function (op::ChebyshevSeriesIntegrationPlan{TF,TI})(
    coeffs::AbstractVector{TF}
) where {TF<:AbstractFloat,TI<:Integer}
    @argcheck length(coeffs) == op.n "length(coeffs) must be equal to n"

    n = length(coeffs)
    tmp = op.tmp
    result = op.result
    weights = op.weights

    # Copy and pad input coefficients
    @inbounds begin
        tmp[1:n] .= coeffs
        tmp[n + 1] = 0
        tmp[n + 2] = 0
    end

    # Compute interior coefficients
    @inbounds begin
        # b₂ = c₁ - c₃/2
        result[2] = tmp[1] - tmp[3] / 2

        # bᵣ = (cᵣ₋₁ - cᵣ₊₁)/(2r) for r > 1
        for r in 2:n
            result[r + 1] = (tmp[r] - tmp[r + 2]) / (2 * r)
        end
    end

    # Compute b₀ to ensure coeffs(-1) = 0
    @inbounds begin
        result[1] = 0
        for i in 1:n
            result[1] += weights[i] * result[i + 1]
        end
    end

    return result
end

function chebyshevt_integrate_plan(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    return ChebyshevSeriesIntegrationPlan{TF}(n)
end

function cheb_series_integrate(coeffs::AbstractVector{TF}) where {TF<:AbstractFloat}
    n = length(coeffs)
    op = ChebyshevSeriesIntegrationPlan{TF}(n)
    return op(coeffs)
end

export cheb_series_integrate
