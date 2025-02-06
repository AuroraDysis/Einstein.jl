"""
    UniformGrid{TF} <: AbstractGrid{TF}

Uniform grid with constant spacing.
"""
struct UniformGrid{TF<:AbstractFloat} <: AbstractGrid{TF}
    x_min::TF # lower bound of the grid
    x_max::TF # upper bound of the grid
    Δx::TF # Grid spacing
    n::Integer # Number of grid points
    grid::StepRangeLen{TF}

    function UniformGrid(x_min::TF, x_max::TF, n::Integer) where {TF<:AbstractFloat}
        @argcheck n >= 0 "n must be nonnegative"
        @argcheck x_max > x_min "x_max must be greater than x_min"

        grid = range(x_min; stop=x_max, length=n)
        dx = step(grid)
        n = length(grid)

        return new{TF}(x_min, x_max, dx, n, grid)
    end
end

Base.length(grid::UniformGrid) = grid.n
Base.step(grid::UniformGrid) = grid.Δx
Base.size(grid::UniformGrid) = (grid.n,)

Base.@propagate_inbounds Base.getindex(grid::UniformGrid, i::Integer) = grid.grid[i]

export UniformGrid

"""
    fdm_grid(x_min::T, x_max::T, dx::T) where {T<:AbstractFloat}
    fdm_grid!(x_grid::Vector{T}, x_min::T, x_max::T, dx::T) where {T<:AbstractFloat}

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

function fdm_grid!(x_grid::Vector{T}, x_min::T, x_max::T, dx::T) where {T<:AbstractFloat}
    @argcheck x_max > x_min "Invalid interval"
    @argcheck dx > 0 "Spacing must be positive"

    n = round(Int, (x_max - x_min) / dx) + 1
    @argcheck length(x_grid) == n "Grid vector length mismatch"

    x_grid_end = x_min + (n - 1) * dx
    @argcheck (x_max - x_grid_end) < 10 * eps(T) "Grid endpoint mismatch: |x_max - x_grid_end| = $(abs(x_max - x_grid_end)) exceeds tolerance ($(10 * eps(T))). Consider adjusting dx to ensure x_max is reached precisely."

    @inbounds for i in 1:n
        x_grid[i] = x_min + (i - 1) * dx
    end

    return x_grid
end

"""
    fdm_check_grid(x_min::T, x_max::T, dx::T) where {T<:AbstractFloat}

Check if the grid spacing is consistent with the interval bounds.

# Arguments
- `x_min`: Lower bound of the interval
- `x_max`: Upper bound of the interval
- `dx`: Grid spacing
"""
function fdm_check_grid(x_min::T, x_max::T, dx::T) where {T<:AbstractFloat}
    n = round(Int, (x_max - x_min) / dx) + 1
    x_grid_end = x_min + (n - 1) * dx
    @argcheck (x_max - x_grid_end) < 10 * eps(T) "Grid endpoint mismatch: |x_max - x_grid_end| = $(abs(x_max - x_grid_end)) exceeds tolerance ($(10 * eps(T))). Consider adjusting dx to ensure x_max is reached precisely."
end

export fdm_grid, fdm_grid!, fdm_check_grid
