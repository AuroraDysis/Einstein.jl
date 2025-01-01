"""
    cheb1_pts([TR=Float64], n::Integer)
    cheb1_pts([TR=Float64], n::Integer, x_min::TR, x_max::TR)

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
x = cheb1_pts(Float64, 5)
x = cheb1_pts(5)  # Same as above

# Generate 5 points mapped to [0,1]
x = cheb1_pts(Float64, 5, 0.0, 1.0)
x = cheb1_pts(5, 0.0, 1.0)  # Same as above
```

See also: [`cheb1_angles`](@ref), [`cheb2_pts`](@ref)
"""
function cheb1_pts(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    @argcheck n >= 0 "n must be nonnegative"

    if n == 0
        return TR[]
    elseif n == 1
        return [zero(TR)]
    end

    x_grid = Vector{TR}(undef, n)

    # Use symmetric indexing for better numerical properties
    pi_over_2n = convert(TR, π) / (2 * n)
    @inbounds begin
        for i in 1:n
            k = -n + 2i - 1
            x_grid[i] = sin(k * pi_over_2n)
        end
    end

    return x_grid
end

function cheb1_pts(n::TI) where {TI<:Integer}
    return cheb1_pts(Float64, n)
end

# Mapped version documentation is inherited from the main docstring
function cheb1_pts(
    ::Type{TR}, n::TI, x_min::TR, x_max::TR
) where {TR<:AbstractFloat,TI<:Integer}
    x_grid = cheb1_pts(TR, n)

    a = (x_max + x_min) / 2
    b = (x_max - x_min) / 2
    @.. x_grid = a + b * x_grid

    return x_grid
end

function cheb1_pts(n::TI, x_min::Float64, x_max::Float64) where {TI<:Integer}
    return cheb1_pts(Float64, n, x_min, x_max)
end

@testset "cheb1_pts" begin
    @testset "n = 5" begin
        n = 5
        x = cheb1_pts(n)

        @test length(x) == n
        @test x ≈ [
            -0.951056516295154,
            -0.587785252292473,
            0.0,
            0.587785252292473,
            0.951056516295154,
        ]
    end

    @testset "n = 6" begin
        n = 6
        x = cheb1_pts(n)

        @test length(x) == n
        @test x ≈ [
            -0.965925826289068,
            -0.707106781186548,
            -0.258819045102521,
            0.258819045102521,
            0.707106781186548,
            0.965925826289068,
        ]
    end
end
