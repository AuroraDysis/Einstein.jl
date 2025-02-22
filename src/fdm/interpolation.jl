"""
    LocalBarycentricInterpolation(points, values; degree=4)

Construct a local barycentric Lagrange interpolant for equispaced data points.

Creates a polynomial approximation of degree `degree` using the data `values` at uniformly 
spaced points `points`. The interpolation is performed locally using `degree + 1` points 
nearest to the evaluation point.

# Arguments
- `points::StepRangeLen{TF}`: Equispaced points for interpolation
- `values::AbstractVector{TF}`: Function values at the points
- `degree::Integer=4`: Degree of the local polynomial interpolant

# Returns
- An interpolant function that can be evaluated at any point within `[minimum(points), maximum(points)]`

# Notes
- Requires `length(points) ≤ degree + 1`
- `points` and `values` must have the same length
- Uses barycentric Lagrange interpolation for numerical stability

# References
- [Berrut2004](@citet*)
"""
struct LocalBarycentricInterpolation{TF<:AbstractFloat}
    points::StepRangeLen{TF}
    values::AbstractVector{TF}
    weights::Vector{TF}
    degree::Integer
    dx_inv::TF

    function LocalBarycentricInterpolation(
        points::StepRangeLen{TF}, values::AbstractVector{TF}; degree::Integer=4
    ) where {TF<:AbstractFloat}
        @argcheck length(points) >= degree + 1 "points is too small for degree"
        @argcheck length(points) == length(values) "points and values must have the same length"

        weights = barycentric_weights(TF, degree)
        dx_inv = inv(step(points))
        return new{TF}(points, values, weights, degree, dx_inv)
    end
end

function (itp::LocalBarycentricInterpolation{TF})(x::TF) where {TF<:AbstractFloat}
    (; points, values, weights, degree, dx_inv) = itp
    @argcheck points[begin] <= x <= points[end] "x is out of range"

    n = length(points)

    # Find indices of nodes nearest to the evaluation point for local interpolation
    relative_pos = (x - points[begin]) * dx_inv  # relative position in points units
    nearest_idx = round(Int, relative_pos)   # index of nearest points point (zero-based)
    num_weights = degree + 1                 # number of points for local interpolation
    half_num_weights, is_odd = divrem(num_weights, 2)

    # Adjust starting index to get the lowest of nearest num_weights points points (one-based)
    start_idx =
        nearest_idx -
        (relative_pos > nearest_idx ? half_num_weights - 1 : half_num_weights - 2 + is_odd)

    # Ensure the interpolation window stays within points boundaries
    start_idx += min(0, n + 1 - start_idx - num_weights) + max(1 - start_idx, 0)
    idx = start_idx:(start_idx + num_weights - 1)
    @inbounds local_points = view(points, idx)
    @inbounds local_values = view(values, idx)

    return barycentric_interpolate(x, local_points, local_values, weights)
end

export LocalBarycentricInterpolation
