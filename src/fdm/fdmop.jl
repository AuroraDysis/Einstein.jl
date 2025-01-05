struct FDMCentralOp{T<:Real}
    der_order::Integer
    acc_order::Integer
    num_coeffs::Integer
    num_side::Integer
    wts::Vector{T}
    one_over_dxn::T
end

# TODO: benchmark this and implement a more efficient version
@inline function (op::FDMCentralOp{T})(
    v::StridedVector{T}
) where {T<:AbstractFloatOrComplex}
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
    wts = fdm_central(T, der_order, acc_order)
    num_coeffs = length(wts)
    num_side = div(num_coeffs - 1, 2)
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
function (op::FDMHermiteOp{T})(
    f::StridedVector{T}, df::StridedVector{T}
) where {T<:AbstractFloatOrComplex}
    return dot(op.Dwts, f) * op.one_over_dxn + dot(op.Ewts, df) * op.one_over_dxnm1
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
    Dwts, Ewts = fdm_hermite(T, der_order, acc_order)
    num_coeffs = length(Dwts)
    num_side = div(num_coeffs - 1, 2)
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

struct FDMDissOp{T<:Real}
    diss_order::Integer
    num_coeffs::Integer
    num_side::Integer
    wts::Vector{T}
    one_over_dx::T
end

# TODO: benchmark this and implement a more efficient version
@inline function (op::FDMDissOp{T})(v::StridedVector{T}) where {T<:AbstractFloatOrComplex}
    return dot(op.wts, v) * op.one_over_dx
end

"""
    fdm_dissop(diss_order::Integer, dx::T) where {T<:AbstractFloat}

Create a finite difference dissipation operator with specified order.

# Arguments
- `diss_order::Integer`: The order of the dissipation operator
- `dx::T`: Grid spacing
"""
function fdm_dissop(diss_order::Integer, dx::T) where {T<:AbstractFloat}
    wts = fdm_disswts(diss_order)
    num_coeffs = length(wts)
    num_side = div(num_coeffs - 1, 2)
    return FDMDissOp{T}(diss_order, num_coeffs, num_side, wts, one(T) / dx)
end

export FDMCentralOp, fdm_centralop
export FDMHermiteOp, fdm_hermiteop
export FDMDissOp, fdm_dissop
