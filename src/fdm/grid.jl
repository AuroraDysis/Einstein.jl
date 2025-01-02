"""
    fdm_grid([TR=Float64], x_min::TR, x_max::TR, dx::TR) where {TR<:AbstractFloat}

Generate a uniform grid for finite difference methods.

# Arguments
- `TR`: Type parameter for the grid points (e.g., Float64)
- `x_min`: Lower bound of the interval
- `x_max`: Upper bound of the interval
- `dx`: Grid spacing

# Returns
- Vector of n uniformly spaced points, where n = round((x_max - x_min)/dx) + 1
"""
function fdm_grid(::Type{TR}, x_min::TR, x_max::TR, dx::TR) where {TR<:AbstractFloat}
    n = round(Int, (x_max - x_min) / dx) + 1

    x_grid_end = x_min + (n - 1) * dx
    @argcheck (x_max - x_grid_end) < 10 * eps(TR) "Grid endpoint mismatch: |x_max - x_grid_end| = $(abs(x_max - x_grid_end)) exceeds tolerance ($(10 * eps(TR))). Consider adjusting dx to ensure x_max is reached precisely."

    x_grid = zeros(TR, n)
    @inbounds for i in 1:n
        x_grid[i] = x_min + (i - 1) * dx
    end

    return x_grid
end

function fdm_grid(x_min::Float64, x_max::Float64, dx::Float64)
    return fdm_grid(Float64, x_min, x_max, dx)
end

export fdm_grid
