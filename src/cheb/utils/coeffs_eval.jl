"""
    cheb_clenshaw(coeffs::AbstractVector{TFC}, x::TFC) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}

Evaluate Chebyshev coefficients at a point using Clenshaw's algorithm.

# Arguments
- `coeffs`: Vector of Chebyshev coefficients ``[c_0, c_1, \\ldots, c_n]``
- `x`: Evaluation point in [-1,1]

# References
- [chebfun/@chebtech/clenshaw.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech/clenshaw.m)
"""
function cheb_clenshaw(
    coeffs::AbstractVector{TFC}, x::TF
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}}
    @argcheck length(coeffs) > 0 "coeffs must have at least one element"

    n = length(coeffs) - 1

    x = 2 * x

    bk1 = zero(TFC)
    bk2 = zero(TFC)

    @inbounds for k in (n + 1):-2:3
        bk2 = coeffs[k] + x * bk1 - bk2
        bk1 = coeffs[k - 1] + x * bk2 - bk1
    end

    # If n is odd, perform the extra step
    if isodd(n)
        tmp = deepcopy(bk1)
        @inbounds bk1 = coeffs[2] + x * bk1 - bk2
        bk2 = tmp
    end

    # Compute the final value
    @inbounds y = coeffs[1] + x * bk1 / 2 - bk2

    return y
end

# TODO: Implement the vectorized version of cheb_clenshaw
function cheb_clenshaw(
    coeffs::AbstractVector{TFC}, x::AbstractArray{TF}
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}}
    return map(xi -> cheb_clenshaw(coeffs, xi), x)
end

"""
    cheb_coeffs_eval(coeffs::AbstractVector{TFC}, x::TF) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}
    cheb_coeffs_eval(coeffs::AbstractVector{TF}, x::AbstractVector{TFC}) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}

Evaluate Chebyshev coefficients at a point.

# Performance Notes
- Clenshaw's algorithm: O(n) operations per point
- (TODO) NDCT: O(n log n) operations for many points simultaneously

# References
- [chebfun/@chebtech/feval.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech/feval.m)
"""
function cheb_coeffs_eval(
    coeffs::AbstractVector{TFC}, x::TF
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}}
    return cheb_clenshaw(coeffs, x)
end

function cheb_coeffs_eval(
    coeffs::AbstractVector{TFC}, x::AbstractVector{TF}
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}}
    return cheb_clenshaw(coeffs, x)
end

export cheb_clenshaw, cheb_coeffs_eval
