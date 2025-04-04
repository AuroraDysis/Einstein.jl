@doc raw"""
    cheb_lobatto_points([TF=Float64], n::Integer) where {TF<:AbstractFloat}
    cheb_lobatto_points([TF=Float64], n::Integer, lower_bound::TF, upper_bound::TF) where {TF<:AbstractFloat}

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
function cheb_lobatto_points(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
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
        x_grid[i + 1] = sin((-nm1 + 2i) * pi_over_2nm1)
    end

    return x_grid
end

function cheb_lobatto_points(n::Integer)
    return cheb_lobatto_points(Float64, n)
end

function cheb_lobatto_points(
    ::Type{TF}, n::Integer, lower_bound::TF, upper_bound::TF
) where {TF<:AbstractFloat}
    x_grid = cheb_lobatto_points(TF, n)

    a = (upper_bound + lower_bound) / 2
    b = (upper_bound - lower_bound) / 2
    @.. x_grid = a + b * x_grid

    return x_grid
end

function cheb_lobatto_points(
    n::Integer, lower_bound::Float64, upper_bound::Float64
)
    return cheb_lobatto_points(Float64, n, lower_bound, upper_bound)
end

export cheb_lobatto_points
