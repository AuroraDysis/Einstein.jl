
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

export GaussChebyshevLobattoGrid
