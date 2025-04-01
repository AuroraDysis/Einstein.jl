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

        points = chebgrid1_points(TF, n, lower_bound, upper_bound)

        return new{TF}(n, lower_bound, upper_bound, points)
    end
end

include("angles.jl")
include("points.jl")

include("coeffs2vals.jl")
include("vals2coeffs.jl")

include("barycentric_weights.jl")
include("quadrature_weights.jl")

include("differentiation_matrix.jl")
include("integration_matrix.jl")

"""
Module for approximating smooth functions on the interval [-1,1] using
function values at 1st-kind Chebyshev points and coefficients of the
corresponding 1st-kind Chebyshev series expansion.
"""
module GaussChebyshevGridTool
const angles = chebgrid1_angles
const points = chebgrid1_points

const coeffs2vals = chebgrid1_coeffs2vals
const vals2coeffs = chebgrid1_vals2coeffs

const barycentric_weights = chebgrid1_barycentric_weights
const quadrature_weights = chebgrid1_quadrature_weights

const differentiation_matrix = chebgrid1_differentiation_matrix
const integration_matrix = chebgrid1_integration_matrix
end

export GaussChebyshevGrid
