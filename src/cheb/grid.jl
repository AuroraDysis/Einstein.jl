"""
Module containing functions for generating Chebyshev grid points.

Chebyshev points are specific sets of points used in polynomial interpolation and
numerical integration. They are particularly useful because they help minimize
Runge's phenomenon and provide excellent interpolation properties.
"""

"""
    ChebGridKind

Enumeration specifying the type of Chebyshev grid points.

# Values
- `FirstKind`: Chebyshev points of the first kind, roots of Chebyshev polynomials
- `SecondKind`: Chebyshev points of the second kind, extrema of Chebyshev polynomials
"""
@enum ChebGridKind begin
    FirstKind = 1
    SecondKind = 2
end

"""
    cheb_grid_first_kind(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}

Generate Chebyshev points of the first kind on the interval [-1, 1].

# Arguments
- `TR`: Type parameter for the grid points (e.g., Float64)
- `n`: Number of grid points

# Returns
- Vector of n Chebyshev points of the first kind

# Mathematical Formula
``x_k = -\\cos\\left(\\frac{(2k + 1)\\pi}{2n}\\right), \\quad k = 0,1,\\ldots,n-1``
"""
function cheb_grid_first_kind(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    grid = zeros(TR, n)

    pi_over_2n = convert(TR, pi) / (2 * n)
    @inbounds begin
        for k in 0:(n - 1)
            grid[k + 1] = -cos((2 * k + 1) * pi_over_2n)
        end
    end

    return grid
end

"""
    cheb_grid_second_kind(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}

Generate Chebyshev points of the second kind on the interval [-1, 1].

# Arguments
- `TR`: Type parameter for the grid points (e.g., Float64)
- `n`: Number of grid points

# Returns
- Vector of n Chebyshev points of the second kind

# Mathematical Formula
``x_k = -\\cos\\left(\\frac{k\\pi}{n-1}\\right), \\quad k = 0,1,\\ldots,n-1``
"""
function cheb_grid_second_kind(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    grid = zeros(TR, n)

    pi_over_nm1 = convert(TR, pi) / (n - 1)
    @inbounds begin
        grid[1] = -1
        grid[n] = 1
        for k in 1:(n - 2)
            grid[k + 1] = -cos(k * pi_over_nm1)
        end
    end

    return grid
end

"""
    cheb_grid(::Type{TR}, n::TI, kind::ChebGridKind=SecondKind) where {TR<:AbstractFloat,TI<:Integer}

Generate Chebyshev grid points of either the first or second kind.

# Arguments
- `TR`: Type parameter for the grid points (e.g., Float64)
- `n`: Number of grid points
- `kind`: Type of Chebyshev points (FirstKind or SecondKind)

# Returns
- Vector of n Chebyshev points of the specified kind

# Examples
```julia
# Generate 10 Chebyshev points of the second kind using Float64
points = cheb_grid(Float64, 10)

# Generate 5 Chebyshev points of the first kind using Float32
points = cheb_grid(Float32, 5, FirstKind)
```
"""
function cheb_grid(
    ::Type{TR}, n::TI, kind::ChebGridKind=SecondKind
) where {TR<:AbstractFloat,TI<:Integer}
    if kind == FirstKind
        return cheb_grid_first_kind(TR, n)
    elseif kind == SecondKind
        return cheb_grid_second_kind(TR, n)
    else
        throw(ArgumentError("kind must be FirstKind or SecondKind"))
    end
end

export cheb_grid
export ChebGridKind, FirstKind, SecondKind
