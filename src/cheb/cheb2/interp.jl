"""
    cheb2_coeffs2vals(coeffs::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}
    op::Cheb2Coeffs2ValsOp([TR=Float64], n::TI)(coeffs::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR},TI<:Integer}

Convert Chebyshev coefficients to values at Chebyshev points of the 2nd kind.

# Performance Guide
For best performance, especially in loops or repeated calls:
```julia
op = Cheb2Coeffs2ValsOp(Float64, n)
values = op(coeffs)
```

# References
- [chebfun/@chebtech2/coeffs2vals.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/coeffs2vals.m)
"""
struct Cheb2InterpOp{TR<:AbstractFloat}
    nodes::Vector{TR}    # Interpolation nodes
    weights::Vector{TR}  # Barycentric weights
    tmp::Vector{TR}      # Temporary storage

    function Cheb2InterpOp(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
        nodes = cheb2_pts(TR, n)
        weights = cheb2_barywts(TR, n)
        tmp = Vector{TR}(undef, n)
        return new{TR}(nodes, weights, tmp)
    end

    function Cheb2InterpOp(n::TI) where {TI<:Integer}
        return Cheb2InterpOp(Float64, n)
    end
end

function (op::Cheb2InterpOp{TR})(
    values::AbstractVector{TR}, x::TR
) where {TR<:AbstractFloat}
    return bary(op.weights, op.nodes, values, x)
end

function cheb2_interp(values::AbstractVector{TR}, x::TR) where {TR<:AbstractFloat}
    n = length(values)
    op = Cheb2InterpOp(TR, n)
    return op(values, x)
end

export cheb2_interp, Cheb2InterpOp

@testset "cheb2_interp" begin
    tol = 100 * eps()

    @testset "Interpolation at nodes" begin
        n = 5
        x = cheb2_pts(n)
        v = sin.(x)
        op = Cheb2InterpOp(n)

        for i in 1:n
            @test op(v, x[i]) ≈ v[i]
        end
    end

    @testset "Interpolation accuracy" begin
        n = 32
        op = Cheb2InterpOp(n)
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
