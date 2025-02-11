@doc raw"""
    cheb1_points([T=Float64], n::Integer) where {T<:AbstractFloat}
    cheb1_points([T=Float64], n::Integer, x_min::T, x_max::T) where {T<:AbstractFloat}

Generate Chebyshev points of the 2nd kind.

For the standard interval [-1,1]:
```math
x_k = -\cos\left(\frac{(2k + 1)\pi}{2n}\right), \quad k = 0,1,\ldots,n-1
```

For mapped interval [x_min,x_max]:
```math
x_{\mathrm{mapped}} = \frac{x_{\mathrm{max}} + x_{\mathrm{min}}}{2} + \frac{x_{\mathrm{min}} - x_{\mathrm{max}}}{2}x_k
```

# Arguments
- `T`: Type parameter for the grid points (e.g., Float64)
- `n`: Number of points
- `x_min`: (Optional) Lower bound of the mapped interval
- `x_max`: (Optional) Upper bound of the mapped interval

# References
- [chebfun/@chebtech1/chebpts.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech1/chebpts.m)
"""
function cheb1_points(::Type{T}, n::Integer) where {T<:AbstractFloat}
    @argcheck n >= 0 "n must be nonnegative"

    if n == 0
        return T[]
    elseif n == 1
        return [zero(T)]
    end

    x_grid = Vector{T}(undef, n)

    # Use symmetric indexing for better numerical properties
    pi_over_2n = convert(T, π) / (2 * n)
    @inbounds begin
        for i in 1:n
            k = -n + 2i - 1
            x_grid[i] = sin(k * pi_over_2n)
        end
    end

    return x_grid
end

function cheb1_points(n::Integer)
    return cheb1_points(Float64, n)
end

# Mapped version
function cheb1_points(::Type{T}, n::Integer, x_min::T, x_max::T) where {T<:AbstractFloat}
    x_grid = cheb1_points(T, n)

    a = (x_max + x_min) / 2
    b = (x_max - x_min) / 2
    @.. x_grid = a + b * x_grid

    return x_grid
end

function cheb1_points(n::Integer, x_min::Float64, x_max::Float64)
    return cheb1_points(Float64, n, x_min, x_max)
end

export cheb1_points
