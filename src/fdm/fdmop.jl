struct FDMCentralOp{T<:Real}
    der_order::Integer
    acc_order::Integer
    num_coeffs::Integer
    num_side::Integer
    wts::Vector{T}
    one_over_dxn::T
end

# TODO: benchmark this and implement a more efficient version
@inline function (op::FDMCentralOp{T})(v::StridedVector{T}) where {T<:AbstractFloat}
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
    return FDMCentralOp{T}(
        der_order, acc_order, num_coeffs, num_side, wts, one(T) / dx^der_order
    )
end

struct FDMHermiteOp{T<:Real}
    der_order::Integer
    acc_order::Integer
    num_coeffs::Integer
    num_side::Integer
    Dwts::Vector{T}
    Ewts::Vector{T}
    one_over_dxn::T
    one_over_dxnm1::T
end

# TODO: benchmark this and implement a more efficient version
function (op::FDMHermiteOp{T})(v::StridedVector{T}) where {T<:AbstractFloat}
    return dot(op.Dwts, v) * op.one_over_dxn + dot(op.Ewts, v) * op.one_over_dxnm1
end

"""
    fdm_hermiteop(der_order::Integer, acc_order::Integer, dx::T) where {T<:AbstractFloat}

Create a Hermite finite difference operator with specified derivative order and accuracy.

# Arguments
- `der_order::Integer`: The order of the derivative to approximate
- `acc_order::Integer`: The desired order of accuracy (must be even)
- `dx::T`: Grid spacing
"""
function fdm_hermiteop(
    der_order::Integer, acc_order::Integer, dx::T
) where {T<:AbstractFloat}
    num_coeffs = fdm_hermitenum(der_order, acc_order)
    num_side = div(num_coeffs - 1, 2)
    Dwts, Ewts = fdm_hermite(T, der_order, acc_order)
    return FDMHermiteOp{T}(
        der_order,
        acc_order,
        num_coeffs,
        num_side,
        Dwts,
        Ewts,
        one(T) / dx^der_order,
        one(T) / dx^(der_order - 1),
    )
end

export FDMCentralOp, fdm_centralop
export FDMHermiteOp, fdm_hermiteop
