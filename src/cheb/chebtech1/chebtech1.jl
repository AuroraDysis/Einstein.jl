include("angles.jl")
include("barycentric_weights.jl")
include("points.jl")
include("quadrature_weights.jl")

include("synthesis.jl")
include("analysis.jl")
include("differentiation_matrix.jl")
include("analysis_matrix.jl")
include("synthesis_matrix.jl")
include("integration_matrix.jl")

module GaussChebyshevGridTool
const angles = chebtech1_angles
const barycentric_weights = chebtech1_barycentric_weights
const points = chebtech1_points
const quadrature_weights = chebtech1_quadrature_weights
end

struct GaussChebyshevGrid{TF<:AbstractFloat} <: AbstractChebyshevGrid{TF}
    n::Integer
    lower_bound::TF
    upper_bound::TF
    points::Vector{TF}

    function GaussChebyshevGrid(
        n::Integer, lower_bound::TF, upper_bound::TF
    ) where {TF<:AbstractFloat}
        @argcheck n >= 0 "n must be nonnegative"
        @argcheck upper_bound > lower_bound "upper_bound must be greater than lower_bound"

        points = chebtech1_points(TF, n, lower_bound, upper_bound)

        return new{TF}(n, lower_bound, upper_bound, points)
    end
end

export GaussChebyshevGrid
