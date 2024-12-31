using FastBroadcast

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

Chebyshev points of the first kind are the roots of Chebyshev polynomials of the first kind.
"""
function cheb_grid_first_kind(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    x_grid = zeros(TR, n)

    pi_over_2n = convert(TR, pi) / (2 * n)
    @inbounds begin
        for k in 0:(n - 1)
            x_grid[k + 1] = -cos((2 * k + 1) * pi_over_2n)
        end
    end

    return x_grid
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

Chebyshev points of the second kind are the extrema of Chebyshev polynomials of the first kind.
"""
function cheb_grid_second_kind(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
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

"""
    cheb_grid(::Type{TR}, n::TI, kind::ChebGridKind=SecondKind) where {TR<:AbstractFloat,TI<:Integer}
    cheb_grid(::Type{TR}, n::TI, x_min::TR, x_max::TR; kind::ChebGridKind=SecondKind) where {TR<:AbstractFloat,TI<:Integer}

Generate Chebyshev grid points of either the first or second kind.

# Arguments
- `TR`: Type parameter for the grid points (e.g., Float64)
- `n`: Number of grid points
- `kind`: Type of Chebyshev points (FirstKind or SecondKind)
- `x_min`: (Optional) Lower bound of the mapped interval
- `x_max`: (Optional) Upper bound of the mapped interval

# Returns
- Vector of n Chebyshev points of the specified kind, either on [-1,1] or mapped to [x_min,x_max]

# Mathematical Details
The standard Chebyshev points are generated on [-1,1]. When mapping to [x_min,x_max],
the following linear transformation is applied:
``x_{mapped} = \\frac{x_{max} + x_{min}}{2} + \\frac{x_{min} - x_{max}}{2}x_{cheb}``

# Examples
```julia
# Generate 10 Chebyshev points of the second kind on [-1,1]
points = cheb_grid(Float64, 10)

# Generate 5 Chebyshev points of the first kind on [-1,1]
points = cheb_grid(Float32, 5, FirstKind)

# Generate 10 Chebyshev points mapped to [0,1]
points = cheb_grid(Float64, 10, 0.0, 1.0)
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

function cheb_grid(
    ::Type{TR}, n::TI, x_min::TR, x_max::TR; kind::ChebGridKind=SecondKind
) where {TR<:AbstractFloat,TI<:Integer}
    x_grid = cheb_grid(TR, n, kind)

    a = (x_max + x_min) / 2
    b = (x_min - x_max) / 2
    @.. x_grid = a + b * x_grid

    return x_grid
end

export cheb_grid
export ChebGridKind, FirstKind, SecondKind
