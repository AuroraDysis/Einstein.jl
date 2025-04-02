@doc raw"""
    cheb_gauss_points([TF=Float64], n::Integer) where {TF<:AbstractFloat}
    cheb_gauss_points([TF=Float64], n::Integer, lower_bound::TF, upper_bound::TF) where {TF<:AbstractFloat}

Generate Chebyshev points of the 2nd kind.

For the standard interval [-1,1]:
```math
x_k = -\cos\left(\frac{(2k + 1)\pi}{2n}\right), \quad k = 0,1,\ldots,n-1
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
- [chebfun/@chebtech1/chebpts.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech1/chebpts.m)
"""
function cheb_gauss_points(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    @argcheck n >= 0 "n must be nonnegative"

    if n == 0
        return TF[]
    elseif n == 1
        return [zero(TF)]
    end

    x_grid = Vector{TF}(undef, n)

    # Use symmetric indexing for better numerical properties
    pi_over_2n = convert(TF, π) / (2 * n)
    @inbounds for i in 1:n
        x_grid[i] = sin((-n + 2i - 1) * pi_over_2n)
    end

    return x_grid
end

function cheb_gauss_points(n::Integer)
    return cheb_gauss_points(Float64, n)
end

function cheb_gauss_points(
    ::Type{TF}, n::Integer, lower_bound::TF, upper_bound::TF
) where {TF<:AbstractFloat}
    x_grid = cheb_gauss_points(TF, n)

    a = (upper_bound + lower_bound) / 2
    b = (upper_bound - lower_bound) / 2
    @.. x_grid = a + b * x_grid

    return x_grid
end

function cheb_gauss_points(n::Integer, lower_bound::Float64, upper_bound::Float64)
    return cheb_gauss_points(Float64, n, lower_bound, upper_bound)
end

export cheb_gauss_points
