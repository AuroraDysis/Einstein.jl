@doc raw"""
    cheb2_points([TF=Float64], n::Integer) where {TF<:AbstractFloat}
    cheb2_points([TF=Float64], n::Integer, x_min::TF, x_max::TF) where {TF<:AbstractFloat}

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
- `TF`: Type parameter for the grid points (e.g., Float64)
- `n`: Number of points
- `x_min`: (Optional) Lower bound of the mapped interval
- `x_max`: (Optional) Upper bound of the mapped interval

# References
- [chebfun/@chebtech2/chebpts.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/chebpts.m)
"""
function cheb2_points(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
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

function cheb2_points(n::Integer)
    return cheb2_points(Float64, n)
end

# Mapped version documentation is inherited from the main docstring
function cheb2_points(
    ::Type{TF}, n::Integer, x_min::TF, x_max::TF
) where {TF<:AbstractFloat}
    x_grid = cheb2_points(TF, n)

    a = (x_max + x_min) / 2
    b = (x_max - x_min) / 2
    @.. x_grid = a + b * x_grid

    return x_grid
end

function cheb2_points(n::Integer, x_min::Float64, x_max::Float64)
    return cheb2_points(Float64, n, x_min, x_max)
end

function _cheb_points(
    ::ChebyshevSecondKindNode, ::Type{TF}, n::Integer, x_min::TF, x_max::TF
) where {TF<:AbstractFloat}
    return cheb2_points(TF, n, x_min, x_max)
end

export cheb2_points
