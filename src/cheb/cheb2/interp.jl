"""
    op::Cheb2InterpOp{TR}(values::Vector{TR}, x::TR) -> TR

Interpolate values at Chebyshev points of the second kind using barycentric interpolation.

# Performance Guide
For best performance, especially in loops or repeated calls:
```julia
# Create operator for size n
op = Cheb2InterpOp{Float64}(n)

# Evaluate at multiple points
for x in points
    result = op(values, x)
end
```

# Mathematical Background
Uses the barycentric interpolation formula with simplified weights for
Chebyshev points of the second kind:
```math
p(x) = \\frac{\\sum_{j=0}^n \\frac{(-1)^j \\delta_j f_j}{x-x_j}}{\\sum_{j=0}^n \\frac{(-1)^j \\delta_j}{x-x_j}}
```
where ``\\delta_j`` is 1/2 at endpoints and 1 elsewhere.

See also: [`cheb2_pts`](@ref), [`cheb2_barywts`](@ref)
"""
struct Cheb2InterpOp{TR<:AbstractFloat}
    nodes::Vector{TR}    # Interpolation nodes
    weights::Vector{TR}  # Barycentric weights
    tmp::Vector{TR}      # Temporary storage

    function Cheb2InterpOp{TR}(n::Integer) where {TR<:AbstractFloat}
        nodes = cheb2_pts(TR, n)
        weights = cheb2_barywts(TR, n)
        tmp = Vector{TR}(undef, n)
        return new{TR}(nodes, weights, tmp)
    end
end

function (op::Cheb2InterpOp{TR})(
    values::AbstractVector{TR}, x::TR
) where {TR<:AbstractFloat}
    return bary(op.weights, op.nodes, values, x)
end

function cheb2_interp(values::AbstractVector{TR}, x::TR) where {TR<:AbstractFloat}
    n = length(values)
    op = Cheb2InterpOp{TR}(n)
    return op(values, x)
end

export cheb2_interp, Cheb2InterpOp

@testset "cheb2_interp" begin
    tol = 100 * eps()

    @testset "Interpolation at nodes" begin
        n = 5
        x = cheb2_pts(n)
        v = sin.(x)
        op = Cheb2InterpOp{Float64}(n)

        for i in 1:n
            @test op(v, x[i]) ≈ v[i]
        end
    end

    @testset "Interpolation accuracy" begin
        n = 32
        op = Cheb2InterpOp{Float64}(n)
        x = cheb2_pts(n)
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
