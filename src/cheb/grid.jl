using FastBroadcast

"""
    cheb1_grid(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    cheb1_grid(::Type{TR}, n::TI, x_min::TR, x_max::TR) where {TR<:AbstractFloat,TI<:Integer}

Generate Chebyshev points of the first kind.

# Arguments
- `TR`: Type parameter for the grid points (e.g., Float64)
- `n`: Number of grid points
- `x_min`: (Optional) Lower bound of the mapped interval
- `x_max`: (Optional) Upper bound of the mapped interval

# Returns
- Vector of n Chebyshev points of the first kind

# Mathematical Details
For the standard interval [-1,1]:
``x_k = -\\cos\\left(\\frac{(2k + 1)\\pi}{2n}\\right), \\quad k = 0,1,\\ldots,n-1``

For mapped interval [x_min,x_max]:
``x_{mapped} = \\frac{x_{max} + x_{min}}{2} + \\frac{x_{min} - x_{max}}{2}x_k``

Chebyshev points of the first kind are the roots of Chebyshev polynomials of the first kind.

# Examples
```julia
# Generate 5 points on [-1,1]
x = cheb1_grid(Float64, 5)

# Generate 5 points mapped to [0,1]
x = cheb1_grid(Float64, 5, 0.0, 1.0)
```
"""
function cheb1_grid(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    x_grid = zeros(TR, n)

    pi_over_2n = convert(TR, pi) / (2 * n)
    @inbounds begin
        for k in 0:(n - 1)
            x_grid[k + 1] = -cos((2 * k + 1) * pi_over_2n)
        end
    end

    return x_grid
end

# Mapped version documentation is inherited from the main docstring
function cheb1_grid(
    ::Type{TR}, n::TI, x_min::TR, x_max::TR
) where {TR<:AbstractFloat,TI<:Integer}
    x_grid = cheb1_grid(TR, n)

    a = (x_max + x_min) / 2
    b = (x_min - x_max) / 2
    @.. x_grid = a + b * x_grid

    return x_grid
end

"""
    cheb2_grid(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    cheb2_grid(::Type{TR}, n::TI, x_min::TR, x_max::TR) where {TR<:AbstractFloat,TI<:Integer}

Generate Chebyshev points of the second kind.

# Arguments
- `TR`: Type parameter for the grid points (e.g., Float64)
- `n`: Number of grid points
- `x_min`: (Optional) Lower bound of the mapped interval
- `x_max`: (Optional) Upper bound of the mapped interval

# Returns
- Vector of n Chebyshev points of the second kind

# Mathematical Details
For the standard interval [-1,1]:
``x_k = -\\cos\\left(\\frac{k\\pi}{n-1}\\right), \\quad k = 0,1,\\ldots,n-1``

For mapped interval [x_min,x_max]:
``x_{mapped} = \\frac{x_{max} + x_{min}}{2} + \\frac{x_{min} - x_{max}}{2}x_k``

Chebyshev points of the second kind are the extrema of Chebyshev polynomials of the first kind.
This includes the endpoints of the interval, making them suitable for boundary value problems.

# Examples
```julia
# Generate 5 points on [-1,1]
x = cheb2_grid(Float64, 5)

# Generate 5 points mapped to [0,π]
x = cheb2_grid(Float64, 5, 0.0, π)
```
"""
function cheb2_grid(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    x_grid = zeros(TR, n)

    pi_over_nm1 = convert(TR, pi) / (n - 1)
    @inbounds begin
        x_grid[1] = -1
        x_grid[n] = 1
        for k in 1:(n - 2)
            x_grid[k + 1] = -cos(k * pi_over_nm1)
        end
    end

    return x_grid
end

# Mapped version documentation is inherited from the main docstring
function cheb2_grid(
    ::Type{TR}, n::TI, x_min::TR, x_max::TR
) where {TR<:AbstractFloat,TI<:Integer}
    x_grid = cheb2_grid(TR, n)

    a = (x_max + x_min) / 2
    b = (x_min - x_max) / 2
    @.. x_grid = a + b * x_grid

    return x_grid
end

export cheb1_grid, cheb2_grid
