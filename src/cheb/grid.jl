abstract type AbstractChebyshevGrid{TF<:AbstractFloat} end

Base.length(grid::AbstractChebyshevGrid) = length(grid.data)
Base.size(grid::AbstractChebyshevGrid) = size(grid.data)

Base.@propagate_inbounds function Base.getindex(grid::AbstractChebyshevGrid, i)
    return grid.data[i]
end

Base.@propagate_inbounds function Base.iterate(grid::AbstractChebyshevGrid, state...)
    return iterate(grid.data, state...)
end

Base.keys(grid::AbstractChebyshevGrid) = keys(grid.data)
Base.firstindex(grid::AbstractChebyshevGrid) = firstindex(grid.data)
Base.lastindex(grid::AbstractChebyshevGrid) = lastindex(grid.data)
Base.eltype(::Type{AbstractChebyshevGrid{TF}}) where {TF} = TF
Base.IndexStyle(::Type{AbstractChebyshevGrid}) = IndexLinear()

struct GaussChebyshevGrid{TF<:AbstractFloat} <: AbstractChebyshevGrid{TF}
    lower_bound::TF
    upper_bound::TF
    data::Vector{TF}

    function GaussChebyshevGrid(
        n::Integer, lower_bound::TF, upper_bound::TF
    ) where {TF<:AbstractFloat}
        @argcheck n >= 0 "n must be nonnegative"
        @argcheck upper_bound > lower_bound "upper_bound must be greater than lower_bound"

        return new{TF}(n, lower_bound, upper_bound)
    end
end

struct GaussChebyshevLobattoGrid{TF<:AbstractFloat} <: AbstractChebyshevGrid{TF}
    lower_bound::TF
    upper_bound::TF
    data::Vector{TF}

    function GaussChebyshevLobattoGrid(
        n::Integer, lower_bound::TF, upper_bound::TF
    ) where {TF<:AbstractFloat}
        @argcheck n >= 0 "n must be nonnegative"
        @argcheck upper_bound > lower_bound "upper_bound must be greater than lower_bound"

        return new{TF}(n, lower_bound, upper_bound)
    end
end

export GaussChebyshevGrid, GaussChebyshevLobattoGrid
