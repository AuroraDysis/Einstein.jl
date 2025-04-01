@doc raw"""
    chebtech2_points([TF=Float64], n::Integer) where {TF<:AbstractFloat}
    chebtech2_points([TF=Float64], n::Integer, lower_bound::TF, upper_bound::TF) where {TF<:AbstractFloat}

Generate Chebyshev points of the 1st kind.

For the standard interval [-1,1]:
```math
x_k = -\cos\left(\frac{k\pi}{n-1}\right), \quad k = 0,1,\ldots,n-1
```

For mapped interval [lower_bound,upper_bound]:
```math
x_{\mathrm{mapped}} = \frac{x_{\mathrm{max}} + x_{\mathrm{min}}}{2} + \frac{x_{\mathrm{min}} - x_{\mathrm{max}}}{2}x_k
```

# Arguments
- `TF`: Type parameter for the grid points (e.g., Float64)
- `n`: Number of points
- `lower_bound`: (Optional) Lower bound of the mapped interval
- `upper_bound`: (Optional) Upper bound of the mapped interval

# References
- [chebfun/@chebtech2/chebpts.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/chebpts.m)
"""
function chebtech2_points(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    @argcheck n >= 0 "n must be nonnegative"

    if n == 0
        return TF[]
    elseif n == 1
        return [zero(TF)]
    end

    x_grid = Vector{TF}(undef, n)

    nm1 = n - 1
    pi_over_2nm1 = convert(TF, π) / (2 * nm1)

    @inbounds for i in 0:nm1
        k = -nm1 + 2i  # Creates range -m:2:m
        x_grid[i + 1] = sin(k * pi_over_2nm1)
    end

    return x_grid
end

function chebtech2_points(n::Integer)
    return chebtech2_points(Float64, n)
end

# Mapped version documentation is inherited from the main docstring
function chebtech2_points(
    ::Type{TF}, n::Integer, lower_bound::TF, upper_bound::TF
) where {TF<:AbstractFloat}
    x_grid = chebtech2_points(TF, n)

    a = (upper_bound + lower_bound) / 2
    b = (upper_bound - lower_bound) / 2
    @.. x_grid = a + b * x_grid

    return x_grid
end

function chebtech2_points(n::Integer, lower_bound::Float64, upper_bound::Float64)
    return chebtech2_points(Float64, n, lower_bound, upper_bound)
end

function _cheb_points(
    ::ChebyshevU, ::Type{TF}, n::Integer, lower_bound::TF, upper_bound::TF
) where {TF<:AbstractFloat}
    return chebtech2_points(TF, n, lower_bound, upper_bound)
end

export chebtech2_points
