"""
    cheb_disswts([TR=Float64], n::TI, a::TI, p::TI) where {TR<:AbstractFloat,TI<:Integer}

Compute exponential dissipation weights for Chebyshev spectral methods [Szilagyi:2009qz](@cite).
```math
w_k = e^{-\\alpha\\left(k / n\\right)^{2 p}}, \\quad k = 0, \\ldots, n-1
```

# Arguments
- `TR`: Type parameter for floating-point precision
- `n`: Number of grid points
- `α`: Dissipation strength parameter
- `p`: Order of dissipation

# Suggested values
- \$\\alpha = 36\$ and \$p = 32\$ for weak dissipation [Szilagyi:2009qz, Hilditch:2015aba](@cite)
- \$\\alpha = 36\$ and \$p = 8\$ for strong dissipation [justin_ripley_2023_8215577](@cite)
"""
function cheb_disswts(::Type{TR}, n::TI, α::TI, p::TI) where {TR<:AbstractFloat,TI<:Integer}
    wts = zeros(TR, n)
    p2 = 2 * p
    ninv = one(TR) / n
    @inbounds for k in 0:(n - 1)
        wts[k + 1] = exp(-α * (k * ninv)^p2)
    end
    return wts
end

function cheb_disswts(n::TI, α::TI, p::TI) where {TI<:Integer}
    return cheb_disswts(Float64, n, α, p)
end

"""
    cheb_dissmat(wts::VT, S::MR1, A::MR2; negsum::Bool = true) where {TR<:AbstractFloat,VT<:AbstractVector{TR},MR1<:AbstractMatrix{TR},MR2<:AbstractMatrix{TR}}

Construct a dissipation matrix using precomputed weights and operators,
optionally applying the 'negative sum trick', which seems make the simulation more stable
according to my tests.

# Arguments
- `wts`: Vector of dissipation weights
- `S`: First operator matrix
- `A`: Second operator matrix
- `negsum`: Boolean flag for negative sum trick (default: true)

The function applies the weights through diagonal scaling and implements
the negative sum trick for diagonal elements when negsum is true.
"""
function cheb_dissmat(
    wts::VT, S::MR1, A::MR2; negsum::Bool=true
) where {
    TR<:AbstractFloat,VT<:AbstractVector{TR},MR1<:AbstractMatrix{TR},MR2<:AbstractMatrix{TR}
}
    F = S * Diagonal(wts) * A

    if !negsum
        return F
    end

    # analog to negative sum trick
    F[diagind(F)] .= 0
    return F[diagind(F)] .+= one(TR) .- sum(F; dims=2)
end
