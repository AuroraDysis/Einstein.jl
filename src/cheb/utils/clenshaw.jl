"""
    cheb_clenshaw(c::AbstractVector{T}, x::T) where {T<:AbstractFloat}

Evaluate Chebyshev coefficients at a point using Clenshaw's algorithm.

# Arguments
- `c`: Vector of Chebyshev coefficients ``[c_0, c_1, \\ldots, c_n]``
- `x`: Evaluation point in [-1,1]

# References
- [chebfun/@chebtech/clenshaw.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech/clenshaw.m)
"""
function cheb_clenshaw(c::AbstractVector{T}, x::T) where {T<:AbstractFloat}
    @argcheck length(c) > 0 "c must have at least one element"

    n = length(c) - 1

    x = 2 * x

    bk1 = zero(T)
    bk2 = zero(T)

    @inbounds for k in (n + 1):-2:3
        bk2 = c[k] + x * bk1 - bk2
        bk1 = c[k - 1] + x * bk2 - bk1
    end

    # If n is odd, perform the extra step
    if isodd(n)
        tmp = deepcopy(bk1)
        @inbounds bk1 = c[2] + x * bk1 - bk2
        bk2 = tmp
    end

    # Compute the final value
    @inbounds y = c[1] + x * bk1 / 2 - bk2

    return y
end

# TODO: Implement the vectorized version of cheb_clenshaw
function cheb_clenshaw(
    c::AbstractVector{TR}, x::AbstractVector{TR}
) where {TR<:AbstractFloat}
    return @inbounds [cheb_clenshaw(c, x[i]) for i in eachindex(x)]
end

export cheb_clenshaw
