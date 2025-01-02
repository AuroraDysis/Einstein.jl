"""
    cheb2_pts([TR=Float64], n::TI) where {TR<:AbstractFloat,TI<:Integer}
    cheb2_pts([TR=Float64], n::TI, x_min::TR, x_max::TR) where {TR<:AbstractFloat,TI<:Integer}

Generate Chebyshev points of the 1st kind.

For the standard interval [-1,1]:
``x_k = -\\cos\\left(\\frac{k\\pi}{n-1}\\right), \\quad k = 0,1,\\ldots,n-1``

For mapped interval [x_min,x_max]:
``x_{mapped} = \\frac{x_{max} + x_{min}}{2} + \\frac{x_{min} - x_{max}}{2}x_k``

# Arguments
- `TR`: Type parameter for the grid points (e.g., Float64)
- `n`: Number of points
- `x_min`: (Optional) Lower bound of the mapped interval
- `x_max`: (Optional) Upper bound of the mapped interval

# References
- [chebfun/@chebtech2/chebpts.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/chebpts.m)
"""
function cheb2_pts(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    @argcheck n >= 0 "n must be nonnegative"

    if n == 0
        return TR[]
    elseif n == 1
        return [zero(TR)]
    end

    x_grid = Vector{TR}(undef, n)

    nm1 = n - 1
    pi_over_2nm1 = convert(TR, π) / (2 * nm1)

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
function cheb2_pts(
    ::Type{TR}, n::TI, x_min::TR, x_max::TR
) where {TR<:AbstractFloat,TI<:Integer}
    x_grid = cheb2_pts(TR, n)

    a = (x_max + x_min) / 2
    b = (x_max - x_min) / 2
    @.. x_grid = a + b * x_grid

    return x_grid
end

function cheb2_pts(n::TI, x_min::Float64, x_max::Float64) where {TI<:Integer}
    return cheb2_pts(Float64, n, x_min, x_max)
end

export cheb2_pts

@testset "cheb2_pts" begin
    @testset "n = 5" begin
        n = 5
        x = cheb2_pts(n)

        @test length(x) == n
        @test x ≈ [-1.0, -0.707106781186548, 0.0, 0.707106781186548, 1.0]
    end

    @testset "n = 6" begin
        n = 6
        x = cheb2_pts(n)

        @test length(x) == n
        @test x ≈ [
            -1.0,
            -0.809016994374948,
            -0.309016994374947,
            0.309016994374947,
            0.809016994374948,
            1.0,
        ]
    end
end
