@doc raw"""
    cheb2_pts([T=Float64], n::Integer) where {T<:AbstractFloat}
    cheb2_pts([T=Float64], n::Integer, x_min::T, x_max::T) where {T<:AbstractFloat}

Generate Chebyshev points of the 1st kind.

For the standard interval [-1,1]:
```math
x_k = -\cos\left(\frac{k\pi}{n-1}\right), \quad k = 0,1,\ldots,n-1
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
- [chebfun/@chebtech2/chebpts.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/chebpts.m)
"""
function cheb2_pts(::Type{T}, n::Integer) where {T<:AbstractFloat}
    @argcheck n >= 0 "n must be nonnegative"

    if n == 0
        return T[]
    elseif n == 1
        return [zero(T)]
    end

    x_grid = Vector{T}(undef, n)

    nm1 = n - 1
    pi_over_2nm1 = convert(T, π) / (2 * nm1)

    @inbounds for i in 0:nm1
        k = -nm1 + 2i  # Creates range -m:2:m
        x_grid[i + 1] = sin(k * pi_over_2nm1)
    end

    return x_grid
end

function cheb2_pts(n::TI) where {TI<:Integer}
    return cheb2_pts(Float64, n)
end

# Mapped version documentation is inherited from the main docstring
function cheb2_pts(::Type{T}, n::Integer, x_min::T, x_max::T) where {T<:AbstractFloat}
    x_grid = cheb2_pts(T, n)

    a = (x_max + x_min) / 2
    b = (x_max - x_min) / 2
    @.. x_grid = a + b * x_grid

    return x_grid
end

function cheb2_pts(n::Integer, x_min::T, x_max::T) where {T<:AbstractFloat}
    return cheb2_pts(T, n, x_min, x_max)
end

export cheb2_pts
