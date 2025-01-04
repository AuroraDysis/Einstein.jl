struct FDMOp{T<:Real}
    der_order::Integer
    acc_order::Integer
    num_coeffs::Integer
    num_side::Integer
    wts::Vector{T}
    one_over_dxn::T
end

# TODO: benchmark this and implement a more efficient version
@inline function (op::FDMOp{T})(v::StridedVector{T}) where {T<:AbstractFloat}
    return dot(op.wts, v) * op.one_over_dxn
end

"""
    fdm_centralop(der_order::Integer, acc_order::Integer, dx::T) where {T<:AbstractFloat}

Create a central finite difference operator with specified derivative order and accuracy.

# Arguments
- `der_order::Integer`: The order of the derivative to approximate
- `acc_order::Integer`: The desired order of accuracy (must be even)
- `dx::T`: Grid spacing
```
"""
function fdm_centralop(
    der_order::Integer, acc_order::Integer, dx::T
) where {T<:AbstractFloat}
    num_coeffs = fdm_centralnum(der_order, acc_order)
    num_side = div(num_coeffs - 1, 2)
    wts = fdm_central(T, der_order, acc_order)
    return FDMOp{T}(der_order, acc_order, num_coeffs, num_side, wts, one(T) / dx^der_order)
end

export FDMOp, fdm_centralop
