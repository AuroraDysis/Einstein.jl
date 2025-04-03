@doc raw"""
    cheb_series_integrate(df::AbstractVector{TF}) where {TF<:AbstractFloat}
    cheb_series_integrate!(f::AbstractVector{TF}, df::AbstractVector{TF}) where {TF<:AbstractFloat}

Compute the indefinite integral of a function $f'(x)$ given its Chebyshev series,
with the constant of integration chosen such that $f(-1) = 0$.

# Arguments
- `df`: Vector of Chebyshev series coefficients of $f'$
- `f`: Vector of Chebyshev series coefficients of $f$

# References
- [chebfun/@chebtech/cumsum.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech/cumsum.m)
"""
function cheb_series_integrate!(
    f::AbstractVector{TF}, df::AbstractVector{TF}
) where {TF<:AbstractFloat}
    @boundscheck begin
        @argcheck length(f) == length(df) + 1 "length(f) must be equal to length(df) + 1"
        @argcheck firstindex(f) == 1 "firstindex(f) must be 1"
        @argcheck firstindex(df) == 1 "firstindex(df) must be 1"
    end

    n = length(df)

    f[2] = df[1] - df[3] / 2
    for r in 2:(n - 2)
        f[r + 1] = (df[r] - df[r + 2]) / (2 * r)
    end
    f[n] = df[n - 1] / (2 * (n - 1))
    f[n + 1] = df[n] / (2 * n)

    b0 = 0
    for r in 1:n
        b0 += (-1)^(r + 1) * f[r + 1]
    end
    f[1] = b0

    return f
end

function cheb_series_integrate(df::AbstractVector{TF}) where {TF<:AbstractFloat}
    f = Array{TF}(undef, length(df) + 1)
    cheb_series_integrate!(f, df)
    return f
end

export cheb_series_integrate, cheb_series_integrate!
