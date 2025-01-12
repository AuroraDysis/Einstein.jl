"""
    fdm_grid(x_min::T, x_max::T, dx::T) where {T<:AbstractFloat}

Generate a uniform grid for finite difference methods.

# Arguments
- `x_min`: Lower bound of the interval
- `x_max`: Upper bound of the interval
- `dx`: Grid spacing

# Returns
- Vector of n uniformly spaced points, where n = round((x_max - x_min)/dx) + 1
"""
function fdm_grid(x_min::T, x_max::T, dx::T) where {T<:AbstractFloat}
    @argcheck x_max > x_min "Invalid interval"
    @argcheck dx > 0 "Spacing must be positive"

    n = round(Int, (x_max - x_min) / dx) + 1

    x_grid_end = x_min + (n - 1) * dx
    @argcheck (x_max - x_grid_end) < 10 * eps(T) "Grid endpoint mismatch: |x_max - x_grid_end| = $(abs(x_max - x_grid_end)) exceeds tolerance ($(10 * eps(T))). Consider adjusting dx to ensure x_max is reached precisely."

    x_grid = Vector{T}(undef, n)
    @inbounds for i in 1:n
        x_grid[i] = x_min + (i - 1) * dx
    end

    return x_grid
end

export fdm_grid
