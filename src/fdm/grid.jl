"""
    UniformGrid{TF} <: AbstractGrid{TF}

Uniform grid with constant spacing.
"""
struct UniformGrid{TF<:AbstractFloat} <: AbstractGrid{TF}
    lower_bound::TF # lower bound of the grid
    upper_bound::TF # upper bound of the grid
    dx::TF # Grid spacing
    data::StepRangeLen{TF}

    function UniformGrid(
        lower_bound::TF, upper_bound::TF, n::Integer
    ) where {TF<:AbstractFloat}
        @argcheck n >= 0 "n must be nonnegative"
        @argcheck upper_bound > lower_bound "upper_bound must be greater than lower_bound"

        data = range(lower_bound; stop=upper_bound, length=n)
        dx = step(data)

        return new{TF}(lower_bound, upper_bound, dx, data)
    end
end

Base.length(grid::UniformGrid) = length(grid.data)
Base.step(grid::UniformGrid) = grid.dx
Base.size(grid::UniformGrid) = size(grid.data)
Base.@propagate_inbounds Base.getindex(grid::UniformGrid, i) = grid.data[i]
Base.keys(grid::UniformGrid) = keys(grid.data)
function Base.iterate(grid::UniformGrid, state=(eachindex(grid.data),))
    return iterate(grid.data, state)
end
Base.firstindex(grid::UniformGrid) = firstindex(grid.data)
Base.lastindex(grid::UniformGrid) = lastindex(grid.data)
Base.eltype(::Type{UniformGrid{TF}}) where {TF<:AbstractFloat} = TF
Base.IndexStyle(::Type{UniformGrid}) = IndexLinear()
Base.diff(grid::UniformGrid) = Fill(grid.dx, length(grid.data) - 1)

"""
    fdm_grid(lower_bound::TF, upper_bound::TF, dx::TF) where {TF<:AbstractFloat}

Create a uniform grid for finite difference methods (FDM).

# Arguments
- `lower_bound::TF`: Lower bound of the grid
- `upper_bound::TF`: Upper bound of the grid
- `dx::TF`: Grid spacing
"""
function fdm_grid(lower_bound::TF, upper_bound::TF, dx::TF) where {TF<:AbstractFloat}
    @argcheck upper_bound > lower_bound "Invalid interval"
    @argcheck dx > 0 "Spacing must be positive"

    n = round(Int, (upper_bound - lower_bound) / dx) + 1
    x_grid_end = lower_bound + (n - 1) * dx
    @argcheck (upper_bound - x_grid_end) < 10 * eps(TF) "Grid endpoint mismatch: |upper_bound - x_grid_end| = $(abs(upper_bound - x_grid_end)) exceeds tolerance ($(10 * eps(TF))). Consider adjusting dx to ensure upper_bound is reached precisely."

    return UniformGrid(lower_bound, upper_bound, n)
end

export fdm_grid, UniformGrid
