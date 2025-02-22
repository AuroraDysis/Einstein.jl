function barycentric_kernal(
    x0::TF, grid::AbstractVector{TF}, values::AbstractVector{TR}, weights::Vector{TF}
) where {TF<:AbstractFloat,TR<:Union{TF,Complex{TF}}}
    p = zero(TR)
    q = zero(TR)

    @inbounds for i in eachindex(grid)
        Δx = x0 - grid[i]

        if iszero(Δx)
            return values[i]
        end

        wi = weights[i] / Δx
        p += wi * values[i]
        q += wi
    end

    return p / q
end

"""
    BarycentricInterpolation{TF<:AbstractFloat}(grid::Vector{TF}, weights::Vector{TF})

A structure representing barycentric interpolation with precomputed weights.

# Fields
- `grid::Vector{TF}`: Vector of interpolation grid (typically Chebyshev grid)
- `weights::Vector{TF}`: Vector of barycentric weights

# Methods
    (itp::BarycentricInterpolation{TF})(y::AbstractVector{TR}, x0::TF) where {TF<:AbstractFloat,TR<:Union{TF,Complex{TF}}}

Evaluate the interpolant at point `x0` for function values `y`.

# Reference
- [chebfun/@chebtech2/bary.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/bary.m)
"""
struct BarycentricInterpolation{TF<:AbstractFloat}
    grid::AbstractVector{TF}   # Grid grid
    weights::Vector{TF}  # Barycentric weights
end

function (itp::BarycentricInterpolation{TF})(
    values::AbstractVector{TR}, x0::TF
) where {TF<:AbstractFloat,TR<:Union{TF,Complex{TF}}}
    (; grid, weights) = itp

    return barycentric_kernal(x0, grid, values, weights)
end

export BarycentricInterpolation
