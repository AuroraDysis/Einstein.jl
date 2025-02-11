"""
    UniformGrid{TF} <: AbstractGrid{TF}

Uniform grid with constant spacing.

# Arguments
- lower_bound::TF: Lower bound of the grid
- upper_bound::TF: Upper bound of the grid
- dx::TF: Grid spacing
"""
struct UniformGrid{TF<:AbstractFloat} <: AbstractGrid{TF}
    lower_bound::TF # lower bound of the grid
    upper_bound::TF # upper bound of the grid
    data::StepRangeLen{TF}

    function UniformGrid(lower_bound::TF, upper_bound::TF, dx::TF) where {TF<:AbstractFloat}
        @argcheck upper_bound > lower_bound "Invalid interval"
        @argcheck dx > 0 "Spacing must be positive"

        n = round(Int, (upper_bound - lower_bound) / dx) + 1
        precise_upper_bound = lower_bound + (n - 1) * dx
        @argcheck (upper_bound - precise_upper_bound) < 10 * eps(TF) "Grid endpoint mismatch: |upper_bound - precise_upper_bound| = $(abs(upper_bound - precise_upper_bound)) exceeds tolerance ($(10 * eps(TF))). Consider adjusting dx to ensure upper_bound is reached precisely."

        data = range(lower_bound, upper_bound, n)

        return new{TF}(lower_bound, upper_bound, data)
    end
end

Base.length(grid::UniformGrid) = length(grid.data)
Base.step(grid::UniformGrid) = step(grid.data)
Base.size(grid::UniformGrid) = size(grid.data)
Base.@propagate_inbounds Base.getindex(grid::UniformGrid, i) = grid.data[i]
Base.keys(grid::UniformGrid) = keys(grid.data)
Base.@propagate_inbounds function Base.iterate(grid::UniformGrid, state...)
    return iterate(grid.data, state...)
end
Base.firstindex(grid::UniformGrid) = firstindex(grid.data)
Base.lastindex(grid::UniformGrid) = lastindex(grid.data)
Base.eltype(::Type{UniformGrid{TF}}) where {TF<:AbstractFloat} = TF
Base.IndexStyle(::Type{UniformGrid}) = IndexLinear()
Base.diff(grid::UniformGrid) = Fill(step(grid.data), length(grid.data) - 1)

export UniformGrid
