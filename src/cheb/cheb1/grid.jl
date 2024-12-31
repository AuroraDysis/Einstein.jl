"""
    cheb1_angle([TR=Float64], n::Integer)

Compute angles for Chebyshev points of the first kind.

# Arguments
- `TR`: Type parameter for the angles (e.g., Float64)
- `n`: Number of points

# Returns
- Vector of n angles in [0,π], ordered decreasing

# Mathematical Details
``\\theta_k = \\frac{(2k + 1)\\pi}{2n}, \\quad k = n-1,\\ldots,0``

These angles generate first-kind Chebyshev points via: x_k = -cos(θ_k)
"""
function cheb1_angle(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    θ = zeros(TR, n)

    pi_over_2n = convert(TR, pi) / (2 * n)
    @inbounds begin
        for k in 0:(n - 1)
            θ[n - k] = (2 * k + 1) * pi_over_2n
        end
    end

    return θ
end

function cheb1_angle(n::TI) where {TI<:Integer}
    return cheb1_angle(Float64, n)
end

"""
    cheb1_grid([TR=Float64], n::Integer)
    cheb1_grid([TR=Float64], n::Integer, x_min::TR, x_max::TR)

Generate Chebyshev points of the first kind.

# Arguments
- `TR`: Type parameter for the grid points (e.g., Float64)
- `n`: Number of points
- `x_min`: (Optional) Lower bound of the mapped interval
- `x_max`: (Optional) Upper bound of the mapped interval

# Returns
- Vector of n Chebyshev points of the first kind

# Mathematical Details
For the standard interval [-1,1]:
``x_k = -\\cos\\left(\\frac{(2k + 1)\\pi}{2n}\\right), \\quad k = 0,1,\\ldots,n-1``

For mapped interval [x_min,x_max]:
``x_{mapped} = \\frac{x_{max} + x_{min}}{2} + \\frac{x_{min} - x_{max}}{2}x_k``

# Notes
Chebyshev points of the first kind are the roots of Chebyshev polynomials T_n(x).
The convenience methods with Integer arguments default to Float64 precision.

# Examples
```julia
# Generate 5 points on [-1,1]
x = cheb1_grid(Float64, 5)
x = cheb1_grid(5)  # Same as above

# Generate 5 points mapped to [0,1]
x = cheb1_grid(Float64, 5, 0.0, 1.0)
x = cheb1_grid(5, 0.0, 1.0)  # Same as above
```

See also: [`cheb1_angle`](@ref), [`cheb2_grid`](@ref)
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

function cheb1_grid(n::TI) where {TI<:Integer}
    return cheb1_grid(Float64, n)
end

# Mapped version documentation is inherited from the main docstring
function cheb1_grid(
    ::Type{TR}, n::TI, x_min::TR, x_max::TR
) where {TR<:AbstractFloat,TI<:Integer}
    x_grid = cheb1_grid(TR, n)

    a = (x_max + x_min) / 2
    b = (x_max - x_min) / 2
    @.. x_grid = a + b * x_grid

    return x_grid
end

function cheb1_grid(n::TI, x_min::Float64, x_max::Float64) where {TI<:Integer}
    return cheb1_grid(Float64, n, x_min, x_max)
end
