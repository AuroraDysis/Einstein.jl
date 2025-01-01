"""
    op::Cheb1InterpOp{TR}(values::Vector{TR}, x::TR) -> TR

Interpolate values at Chebyshev points of the first kind using barycentric interpolation.

# Performance Guide
For best performance, especially in loops or repeated calls:
```julia
# Create operator for size n
op = Cheb1InterpOp{Float64}(n)

# Evaluate at multiple points
for x in points
    result = op(values, x)
end
```

# Mathematical Background
Uses the barycentric interpolation formula:
```math
p(x) = \\frac{\\sum_{j=0}^n \\frac{w_j f_j}{x-x_j}}{\\sum_{j=0}^n \\frac{w_j}{x-x_j}}
```
where ``w_j`` are the barycentric weights and ``x_j`` are the Chebyshev points.

See also: [`cheb1_pts`](@ref), [`cheb1_barywts`](@ref)
"""
struct Cheb1InterpOp{TR<:AbstractFloat}
    nodes::Vector{TR}    # Interpolation nodes
    weights::Vector{TR}  # Barycentric weights
    tmp::Vector{TR}      # Temporary storage

    function Cheb1InterpOp{TR}(n::Integer) where {TR<:AbstractFloat}
        nodes = cheb1_pts(TR, n)
        weights = cheb1_barywts(TR, n)
        tmp = Vector{TR}(undef, n)
        return new{TR}(nodes, weights, tmp)
    end
end

function (op::Cheb1InterpOp{TR})(
    values::AbstractVector{TR}, x::TR
) where {TR<:AbstractFloat}
    return bary(op.weights, op.nodes, values, x)
end

function cheb1_interp(values::AbstractVector{TR}, x::TR) where {TR<:AbstractFloat}
    n = length(values)
    op = Cheb1InterpOp{TR}(n)
    return op(values, x)
end

export cheb1_interp, Cheb1InterpOp

@testset "cheb1_interp" begin
    tol = 100 * eps()

    @testset "Interpolation at nodes" begin
        n = 5
        x = cheb1_pts(n)
        v = sin.(x)
        op = Cheb1InterpOp{Float64}(n)

        for i in 1:n
            @test op(v, x[i]) ≈ v[i]
        end
    end

    @testset "Interpolation accuracy" begin
        n = 32
        op = Cheb1InterpOp{Float64}(n)
        x = cheb1_pts(n)
        v = cos.(x)

        # Test at intermediate points
        test_points = [-0.9, -0.5, 0.0, 0.5, 0.9]
        for xi in test_points
            y_interp = op(v, xi)
            y_exact = cos(xi)
            @test abs(y_interp - y_exact) < tol
        end
    end
end
