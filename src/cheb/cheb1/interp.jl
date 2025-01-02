"""
    cheb1_interp(values::VT, x::TR) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}
    Cheb1InterpOp([TR=Float64], n::TI)(values::VT, x::TR) where {TR<:AbstractFloat,TI<:Integer,VT<:AbstractVector{TR}}

Interpolate values at Chebyshev points of the 1st kind using barycentric interpolation.

# Performance Guide
For best performance, especially in loops or repeated calls:
```julia
op = Cheb1InterpOp(Float64, n)
y = op(v, x)
```
"""
struct Cheb1InterpOp{TR<:AbstractFloat}
    nodes::Vector{TR}    # Interpolation nodes
    weights::Vector{TR}  # Barycentric weights
    tmp::Vector{TR}      # Temporary storage

    function Cheb1InterpOp(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
        nodes = cheb1_pts(TR, n)
        weights = cheb1_barywts(TR, n)
        tmp = Vector{TR}(undef, n)
        return new{TR}(nodes, weights, tmp)
    end

    function Cheb1InterpOp(n::TI) where {TI<:Integer}
        return Cheb1InterpOp(Float64, n)
    end
end

function (op::Cheb1InterpOp{TR})(
    values::AbstractVector{TR}, x::TR
) where {TR<:AbstractFloat}
    return bary(op.weights, op.nodes, values, x)
end

function cheb1_interp(values::AbstractVector{TR}, x::TR) where {TR<:AbstractFloat}
    n = length(values)
    op = Cheb1InterpOp(TR, n)
    return op(values, x)
end

export cheb1_interp, Cheb1InterpOp

@testset "cheb1_interp" begin
    tol = 100 * eps()

    @testset "Interpolation at nodes" begin
        n = 5
        x = cheb1_pts(n)
        v = sin.(x)
        op = Cheb1InterpOp(n)

        for i in 1:n
            @test op(v, x[i]) ≈ v[i]
        end
    end

    @testset "Interpolation accuracy" begin
        n = 32
        op = Cheb1InterpOp(n)
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
