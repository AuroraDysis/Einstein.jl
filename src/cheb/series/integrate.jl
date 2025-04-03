@doc raw"""
    cheb_series_integrate(df::AbstractVector{TFC}) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    cheb_series_integrate!(f::AbstractVector{TFC}, df::AbstractVector{TFC}) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}

Compute the indefinite integral of a function $f'(x)$ given its Chebyshev series,
with the constant of integration chosen such that $f(-1) = 0$.

# Mathematical Details

If the input function $f'(x)$ is represented by a Chebyshev series of length $n$:
```math
f'(x) = \sum_{r=0}^{n-1} c_r T_r(x)
```
its integral $f(x)$ is represented by a Chebyshev series of length $n+1$:
```math
f(x) = \sum_{r=0}^{n} b_r T_r(x)
```
where:
- $b_0$ is determined from the constant of integration as:
```math
b_0 = \sum_{r=1}^{n} (-1)^{r+1} b_r
```
- The other coefficients are given by:
```math
b_1 = c_0 - c_2/2
b_r = (c_{r-1} - c_{r+1})/(2r) \text{ for } r > 1
```
with $c_{n+1} = c_{n+2} = 0$.

# Arguments
- `df`: Vector of Chebyshev series coefficients of $f'$
- `f`: Vector of Chebyshev series coefficients of $f$

# References
- [chebfun/@chebtech/cumsum.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech/cumsum.m)
"""
function cheb_series_integrate!(
    f::AbstractVector{TFC}, df::AbstractVector{TFC}
) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    @boundscheck begin
        @argcheck length(f) == length(df) + 1 "length(f) must be equal to length(df) + 1"
        @argcheck length(df) >= 1 "length(df) must be greater than 1"
        @argcheck firstindex(f) == 1 "firstindex(f) must be 1"
        @argcheck firstindex(df) == 1 "firstindex(df) must be 1"
    end

    n = length(df)

    if n == 1
        f[1] = df[1]
        f[2] = df[1]
        return f
    elseif n == 2
        f[1] = df[1] - df[2] / 4
        f[2] = df[1]
        f[3] = df[2] / 4
        return f
    end

    begin
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
    end

    return f
end

function cheb_series_integrate(
    df::AbstractVector{TFC}
) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    f = Array{TFC}(undef, length(df) + 1)
    cheb_series_integrate!(f, df)
    return f
end

export cheb_series_integrate, cheb_series_integrate!
