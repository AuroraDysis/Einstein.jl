"""
    fdm_grid(::Type{TR}, x_min::TR, x_max::TR, dx::TR) where {TR<:AbstractFloat}

Generate a uniform grid for finite difference methods.

# Arguments
- `TR`: Type parameter for the grid points (e.g., Float64)
- `x_min`: Lower bound of the interval
- `x_max`: Upper bound of the interval
- `dx`: Grid spacing

# Returns
- Vector of n uniformly spaced points, where n = round((x_max - x_min)/dx) + 1

# Mathematical Details
The grid points are generated as:
```math
x_i = x_{min} + (i-1)dx, \\quad i = 1,\\ldots,n
```
where n is chosen to ensure the grid covers [x_min, x_max] with spacing dx.

# Notes
- The function ensures that x_max is accurately represented within floating-point precision
- The actual number of points is computed to maintain uniform spacing
- The final point may differ from x_max by at most machine epsilon

# Examples
```julia
# Generate grid with spacing 0.1 on [0,1]
x = fdm_grid(Float64, 0.0, 1.0, 0.1)

# Generate grid with 100 points on [-1,1]
dx = 2.0/99  # To get exactly 100 points
x = fdm_grid(Float64, -1.0, 1.0, dx)
```
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

export fdm_grid
