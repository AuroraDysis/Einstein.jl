"""
    UniformGrid{TF} <: AbstractGrid{TF}

Uniform grid with constant spacing.
"""
struct UniformGrid{TF<:AbstractFloat} <: AbstractGrid{TF}
    x_min::TF # lower bound of the grid
    x_max::TF # upper bound of the grid
    Δx::TF # Grid spacing
    n::Integer # Number of grid points
    data::StepRangeLen{TF}

    function UniformGrid(x_min::TF, x_max::TF, n::Integer) where {TF<:AbstractFloat}
        @argcheck n >= 0 "n must be nonnegative"
        @argcheck x_max > x_min "x_max must be greater than x_min"

        data = range(x_min; stop=x_max, length=n)
        dx = step(data)
        n = length(data)

        return new{TF}(x_min, x_max, dx, n, data)
    end
end

Base.length(grid::UniformGrid) = grid.n
Base.step(grid::UniformGrid) = grid.Δx
Base.size(grid::UniformGrid) = (grid.n,)
Base.@propagate_inbounds Base.getindex(grid::UniformGrid, i) = grid.data[i]
Base.keys(grid::UniformGrid) = keys(grid.data)
Base.iterate(grid::UniformGrid, state=(eachindex(grid.data),)) = iterate(grid.data, state)
Base.firstindex(grid::UniformGrid) = firstindex(grid.data)
Base.lastindex(grid::UniformGrid) = lastindex(grid.data)
Base.eltype(::Type{UniformGrid{TF}}) where {TF<:AbstractFloat} = TF
Base.IndexStyle(::Type{UniformGrid}) = IndexLinear()
Base.diff(grid::UniformGrid) = Fill(grid.Δx, grid.n - 1)

"""
    fdm_grid(x_min::TF, x_max::TF, dx::TF) where {TF<:AbstractFloat}

Create a uniform grid for finite difference methods (FDM).

# Arguments
- `x_min::TF`: Lower bound of the grid
- `x_max::TF`: Upper bound of the grid
- `dx::TF`: Grid spacing
"""
function fdm_grid(x_min::TF, x_max::TF, dx::TF) where {TF<:AbstractFloat}
    @argcheck x_max > x_min "Invalid interval"
    @argcheck dx > 0 "Spacing must be positive"

    n = round(Int, (x_max - x_min) / dx) + 1
    x_grid_end = x_min + (n - 1) * dx
    @argcheck (x_max - x_grid_end) < 10 * eps(TF) "Grid endpoint mismatch: |x_max - x_grid_end| = $(abs(x_max - x_grid_end)) exceeds tolerance ($(10 * eps(TF))). Consider adjusting dx to ensure x_max is reached precisely."

    return UniformGrid(x_min, x_max, n)
end

export fdm_grid, UniformGrid
