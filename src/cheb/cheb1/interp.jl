"""
    cheb1_interp(values::VT, x::TR) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}
    Cheb1InterpOp([TR=Float64], n::Integer)(values::VT, x::TR) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}

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

    function Cheb1InterpOp(::Type{TR}, n::Integer) where {TR<:AbstractFloat}
        nodes = cheb1_pts(TR, n)
        weights = cheb1_barywts(TR, n)
        tmp = Vector{TR}(undef, n)
        return new{TR}(nodes, weights, tmp)
    end

    function Cheb1InterpOp(n::Integer)
        return Cheb1InterpOp(Float64, n)
    end
end

function (op::Cheb1InterpOp{TR})(
    values::AbstractVector{TR}, x::TR
) where {TR<:AbstractFloat}
    return bary(op.weights, op.nodes, values, x)
end

function cheb1_interp(values::VT, x::TR) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}
    n = length(values)
    op = Cheb1InterpOp(TR, n)
    return op(values, x)
end

export cheb1_interp, Cheb1InterpOp
