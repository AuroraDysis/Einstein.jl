"""
    cheb1_pts([TR=Float64], n::TI) where {TR<:AbstractFloat,TI<:Integer}
    cheb1_pts([TR=Float64], n::TI, x_min::TR, x_max::TR) where {TR<:AbstractFloat,TI<:Integer}

Generate Chebyshev points of the 2nd kind.

For the standard interval [-1,1]:
``x_k = -\\cos\\left(\\frac{(2k + 1)\\pi}{2n}\\right), \\quad k = 0,1,\\ldots,n-1``

For mapped interval [x_min,x_max]:
``x_{mapped} = \\frac{x_{max} + x_{min}}{2} + \\frac{x_{min} - x_{max}}{2}x_k``

# Arguments
- `TR`: Type parameter for the grid points (e.g., Float64)
- `n`: Number of points
- `x_min`: (Optional) Lower bound of the mapped interval
- `x_max`: (Optional) Upper bound of the mapped interval

# References
- [chebfun/@chebtech1/chebpts.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech1/chebpts.m)
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

export cheb1_pts

@testset "cheb1_pts" begin
    for TR in [Float64, BigFloat]
        tol = 10 * eps(TR)
        for intervals in [[-1, 1], [0, 1]]
            x_min = convert(TR, intervals[1])
            x_max = convert(TR, intervals[2])
            @testset "$TR, $intervals" begin
                for n in 0:10
                    @test isapprox(
                        cheb1_pts(TR, n, x_min, x_max),
                        reverse(points(Chebyshev(x_min .. x_max), n)),
                        atol=tol,
                    )
                end
            end
        end
    end
end
