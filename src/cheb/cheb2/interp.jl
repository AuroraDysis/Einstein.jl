"""
    cheb2_coeffs2vals(coeffs::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}
    op::Cheb2Coeffs2ValsOp([TR=Float64], n::Integer)(coeffs::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}

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

    function Cheb2InterpOp(::Type{TR}, n::Integer) where {TR<:AbstractFloat}
        nodes = cheb2_pts(TR, n)
        weights = cheb2_barywts(TR, n)
        tmp = Vector{TR}(undef, n)
        return new{TR}(nodes, weights, tmp)
    end

    function Cheb2InterpOp(n::Integer)
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
