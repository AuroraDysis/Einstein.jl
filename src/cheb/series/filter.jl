@doc raw"""
    cheb_series_filter_weights_exp([T=Float64], n::Integer, a::Integer, p::Integer) where {T<:AbstractFloat}

Compute exponential filter weights for Chebyshev series [Szilagyi:2009qz](@cite).
```math
w_k = e^{-\alpha\left(k / n\right)^{2 p}}, \quad k = 0, \ldots, n-1
```

# Arguments
- `T`: Type parameter for floating-point precision
- `n`: Number of grid points
- `α`: Filter strength parameter
- `p`: Order of filter

# Examples
- $\alpha = 36$ and $p = 32$ for weak filter [Szilagyi:2009qz, Hilditch:2015aba](@cite)
- $\alpha = 40$ and $p = 8$ for strong filter [justin_ripley_2023_8215577](@cite)
"""
function cheb_series_filter_weights_exp(
    ::Type{T}, n::Integer, α::Integer, p::Integer
) where {T<:AbstractFloat}
    wts = zeros(T, n)
    p2 = 2 * p
    ninv = one(T) / n
    @inbounds for k in 0:(n - 1)
        wts[k + 1] = exp(-α * (k * ninv)^p2)
    end
    return wts
end

function cheb_series_filter_weights_exp(n::Integer, α::Integer, p::Integer)
    return cheb_series_filter_weights_exp(Float64, n, α, p)
end

"""
    cheb_filter_matrix(weights::AbstractVector{T}, S::AbstractMatrix{T}, A::AbstractMatrix{T}; negsum::Bool=true) where T<:AbstractFloat

Construct a filter matrix using precomputed weights and operators
for Chebyshev collocation methods, optionally applying the 'negative sum trick',
which seems make the simulation more stable according to my tests.

# Arguments
- `weights`: Vector of filter weights
- `S`: First operator matrix
- `A`: Second operator matrix
- `negative_sum_trick`: Boolean flag for negative sum trick (default: true)

The function applies the weights through diagonal scaling and implements
the negative sum trick for diagonal elements when `negative_sum_trick` is true.
"""
function cheb_filter_matrix(
    weights::AbstractVector{T},
    S::AbstractMatrix{T},
    A::AbstractMatrix{T};
    negative_sum_trick::Bool=true,
) where {T<:AbstractFloat}
    F = S * Diagonal(weights) * A

    if !negative_sum_trick
        return F
    end

    # analog to negative sum trick
    F[diagind(F)] .= 0
    F[diagind(F)] .+= one(T) .- sum(F; dims=2)

    return F
end

export cheb_series_filter_weights_exp, cheb_filter_matrix
