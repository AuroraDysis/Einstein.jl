"""
    cheb_cumsum(f::AbstractVector{TR}) where {TR<:AbstractFloat}
    ChebCumsumOp{[TR=Float64]}(n::TI)(f::AbstractVector{TR}) where {TR<:AbstractFloat,TI<:Integer}

Compute the indefinite integral of a function given its Chebyshev coefficients.

# Arguments
- `f`: Vector of Chebyshev coefficients of the function to be integrated

# References
- [chebfun/@chebtech/cumsum.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech/cumsum.m)
"""
struct ChebCumsumOp{TR<:AbstractFloat,TI<:Integer}
    n::TI              # Number of coefficients
    tmp::Vector{TR}    # Temporary storage for padded coefficients
    result::Vector{TR} # Result storage
    v::Vector{TI}      # Pre-computed alternating signs

    function ChebCumsumOp{TR}(n::TI) where {TR<:AbstractFloat,TI<:Integer}
        # Pre-allocate workspace
        tmp = Vector{TR}(undef, n + 2)
        result = Vector{TR}(undef, n + 1)

        # Pre-compute alternating signs [1, -1, 1, -1, ...]
        v = ones(TI, n)
        @inbounds for i in 2:2:n
            v[i] = -one(TI)
        end

        return new{TR,TI}(n, tmp, result, v)
    end

    function ChebCumsumOp(n::TI) where {TI<:Integer}
        return ChebCumsumOp{Float64}(n)
    end
end

function (op::ChebCumsumOp{TR,TI})(
    f::AbstractVector{TR}
) where {TR<:AbstractFloat,TI<:Integer}
    @argcheck length(f) == op.n "length(f) must be equal to n"

    n = length(f)
    tmp = op.tmp
    result = op.result
    v = op.v

    # Copy and pad input coefficients
    @inbounds begin
        tmp[1:n] .= f
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

    # Compute b₀ to ensure f(-1) = 0
    @inbounds begin
        result[1] = 0
        for i in 1:n
            result[1] += v[i] * result[i + 1]
        end
    end

    return result
end

function cheb_cumsum(f::AbstractVector{TR}) where {TR<:AbstractFloat}
    n = length(f)
    op = ChebCumsumOp{TR}(n)
    return op(f)
end

export cheb_cumsum, ChebCumsumOp
