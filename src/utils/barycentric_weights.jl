@doc raw"""
    barycentric_weights(x::AbstractVector{TF}) where {TF<:AbstractFloat}
    barycentric_weights(x::AbstractRange{TF}) where {TF<:AbstractFloat}
    barycentric_weights([TF=Float64], order::Integer) where {TF<:AbstractFloat}

Compute normalized barycentric weights for interpolation nodes. These weights are used in barycentric Lagrange interpolation.
For Chebyshev-Gauss nodes or Chebyshev-Lobatto nodes,
[cheb_gauss_barycentric_weights](@ref Einstein.ChebyshevSuite.cheb_gauss_barycentric_weights)
and [cheb_lobatto_barycentric_weights](@ref Einstein.ChebyshevSuite.cheb_lobatto_barycentric_weights)
are more efficient implementations.

# Arguments
- `x::AbstractVector{TF}`: Vector of distinct interpolation nodes, for equidistant nodes, use range is recommended `x = x_min:dx:x_max`
- `order::Integer`: order of the interpolation (order = length(x) - 1)

# Returns
- `Vector{TF}`: Barycentric weights scaled such that their infinity norm equals 1

# Details
The barycentric weights wⱼ for nodes xⱼ are computed as:
    wⱼ = 1/∏(xⱼ - xₖ) for k ≠ j

For equidistant nodes, a more efficient formula is used based on binomial coefficients.

# Examples
```julia
# For arbitrary nodes
x = [0.0, 1.0, 2.0]
w = barycentric_weights(x)

# For equidistant nodes
w = barycentric_weights(2)  # 3 nodes: 0, 1, 2

# or
x = 0.0:0.1:1.0
w = barycentric_weights(x)
```

# Notes
- The weights are scaled to have unit infinity norm for numerical stability
- Returns NaN weights if input points are not distinct
- log-sum-exp trick is used to prevent underflow/overflow
- For equidistant nodes, a more efficient algorithm based on binomial coefficients is used

# References
- [Berrut2004](@citet*)
- [chebfun/baryWeights.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/baryWeights.m)
"""
function barycentric_weights(x::AbstractVector{TF}) where {TF<:AbstractFloat}
    @boundscheck begin
        @argcheck length(x) > 1 "x must have at least two elements"
    end

    weights = similar(x)

    x_min, x_max = extrema(x)
    if x_min ≈ x_max
        throw(ArgumentError("Input points are not distinct. Returning NaN weights."))
    end

    # Capacity of the interval - helps prevent underflow/overflow
    c = 4 / (x_max - x_min)

    diff = similar(x)
    for i in eachindex(x)
        @.. diff = c * (x[i] - x)
        diff[i] = 1

        # Log-Sum-Exp trick to avoid underflow/overflow
        prod_sign = prod(sign.(diff))
        @.. diff = log(abs(diff))
        weights[i] = 1 / (prod_sign * exp(sum(diff)))
    end

    weights ./= norm(weights, Inf)

    return weights
end

function barycentric_weights(x::AbstractRange{TF}) where {TF<:AbstractFloat}
    order = length(x) - 1
    return barycentric_weights(TF, order)
end

function barycentric_weights(::Type{TF}, order::Integer) where {TF<:AbstractFloat}
    @boundscheck begin
        @argcheck order > 0 "order must be a positive integer"
    end

    weights = Array{TF}(undef, order + 1)
    half_order = order ÷ 2
    norm_factor, _ = logabsbinomial(order, half_order)
    for j in 0:order
        s = isodd(j) ? -1 : 1
        log_nj, _ = logabsbinomial(order, j)
        weights[j + 1] = s * exp(log_nj - norm_factor)
    end

    return weights
end

function barycentric_weights(order::Integer)
    return barycentric_weights(Float64, order)
end

export barycentric_weights
