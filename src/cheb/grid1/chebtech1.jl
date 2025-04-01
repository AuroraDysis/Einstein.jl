include("angles.jl")
include("points.jl")

include("coeffs2vals.jl")
include("vals2coeffs.jl")

include("barycentric_weights.jl")
include("quadrature_weights.jl")

include("differentiation_matrix.jl")
include("analysis_matrix.jl")
include("synthesis_matrix.jl")
include("integration_matrix.jl")

"""
Module for approximating smooth functions on the interval [-1,1] using
function values at 1st-kind Chebyshev points and coefficients of the
corresponding 1st-kind Chebyshev series expansion.
"""
module GaussChebyshevGridTool
const angles = chebtech1_angles
const points = chebtech1_points

const coeffs2vals = chebtech1_coeffs2vals
const vals2coeffs = chebtech1_vals2coeffs

const barycentric_weights = chebtech1_barycentric_weights
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
