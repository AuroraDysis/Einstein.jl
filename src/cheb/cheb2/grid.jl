"""
    cheb2_angle([TR=Float64], n::Integer)

Compute angles for Chebyshev points of the second kind.

# Arguments
- `TR`: Type parameter for the angles (e.g., Float64)
- `n`: Number of points

# Returns
- Vector of n angles in [0,π], ordered decreasing
- Empty vector if n=0
- [π/2] if n=1

# Mathematical Details
For n > 1:
``\\theta_k = \\frac{k\\pi}{n-1}, \\quad k = n-1,\\ldots,0``

These angles generate second-kind Chebyshev points via: x_k = -cos(θ_k)
"""
function cheb2_angle(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    @argcheck n >= 0 "n must be nonnegative"

    if n == 0
        return TR[]
    elseif n == 1
        return TR[convert(TR, π) / 2]
    end

    θ = zeros(TR, n)
    nm1 = n - 1
    pi_over_nm1 = convert(TR, π) / nm1

    @inbounds for k in 0:nm1
        θ[n - k] = k * pi_over_nm1
    end

    return θ
end

function cheb2_angle(n::TI) where {TI<:Integer}
    return cheb2_angle(Float64, n)
end

"""
    cheb2_grid([TR=Float64], n::Integer)
    cheb2_grid([TR=Float64], n::Integer, x_min::TR, x_max::TR)

Generate Chebyshev points of the second kind.

# Arguments
- `TR`: Type parameter for the grid points (e.g., Float64)
- `n`: Number of points
- `x_min`: (Optional) Lower bound of the mapped interval
- `x_max`: (Optional) Upper bound of the mapped interval

# Returns
- Vector of n Chebyshev points of the second kind

# Mathematical Details
For the standard interval [-1,1]:
``x_k = -\\cos\\left(\\frac{k\\pi}{n-1}\\right), \\quad k = 0,1,\\ldots,n-1``

For mapped interval [x_min,x_max]:
``x_{mapped} = \\frac{x_{max} + x_{min}}{2} + \\frac{x_{min} - x_{max}}{2}x_k``

# Notes
Chebyshev points of the second kind are the extrema of Chebyshev polynomials T_n(x).
These points include the endpoints, making them suitable for boundary value problems.
The convenience methods with Integer arguments default to Float64 precision.

# Examples
```julia
# Generate 5 points on [-1,1]
x = cheb2_grid(Float64, 5)
x = cheb2_grid(5)  # Same as above

# Generate 5 points mapped to [0,π]
x = cheb2_grid(Float64, 5, 0.0, π)
x = cheb2_grid(5, 0.0, π)  # Same as above
```

See also: [`cheb2_angle`](@ref), [`cheb1_grid`](@ref)
"""
function cheb2_grid(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    @argcheck n >= 0 "n must be nonnegative"

    if n == 0
        return TR[]
    elseif n == 1
        return [zero(TR)]
    end

    x_grid = zeros(TR, n)

    nm1 = n - 1
    pi_over_2nm1 = convert(TR, π) / (2 * nm1)

    @inbounds for i in 0:nm1
        k = -nm1 + 2i  # Creates range -m:2:m
        x_grid[i + 1] = sin(k * pi_over_2nm1)
    end

    return x_grid
end

function cheb2_grid(n::TI) where {TI<:Integer}
    return cheb2_grid(Float64, n)
end

# Mapped version documentation is inherited from the main docstring
function cheb2_grid(
    ::Type{TR}, n::TI, x_min::TR, x_max::TR
) where {TR<:AbstractFloat,TI<:Integer}
    x_grid = cheb2_grid(TR, n)

    a = (x_max + x_min) / 2
    b = (x_max - x_min) / 2
    @.. x_grid = a + b * x_grid

    return x_grid
end

function cheb2_grid(n::TI, x_min::Float64, x_max::Float64) where {TI<:Integer}
    return cheb2_grid(Float64, n, x_min, x_max)
end

export cheb2_angle, cheb2_grid
