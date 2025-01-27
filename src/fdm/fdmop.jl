struct FDMCentralOp{TR<:Real}
    der_order::Integer
    acc_order::Integer
    num_coeffs::Integer
    num_side::Integer
    wts::Vector{TR}
    one_over_dxn::TR
end

# TODO: benchmark this and implement a more efficient version
@inline function (op::FDMCentralOp{TR})(
    v::StridedVector{TRC}
) where {TR<:AbstractFloat,TRC<:Union{TR,Complex{TR}}}
    return dot(op.wts, v) * op.one_over_dxn
end

"""
    fdm_centralop(der_order::Integer, acc_order::Integer, dx::TR) where {TR<:AbstractFloat}

Create a central finite difference operator with specified derivative order and accuracy.

# Arguments
- `der_order::Integer`: The order of the derivative to approximate
- `acc_order::Integer`: The desired order of accuracy (must be even)
- `dx::TR`: Grid spacing
```
"""
function fdm_centralop(
    der_order::Integer, acc_order::Integer, dx::TR
) where {TR<:AbstractFloat}
    wts = fdm_centralwts(TR, der_order, acc_order)
    num_coeffs = length(wts)
    num_side = div(num_coeffs - 1, 2)
    return FDMCentralOp{TR}(
        der_order, acc_order, num_coeffs, num_side, wts, one(TR) / dx^der_order
    )
end

struct FDMBoundOp{TR<:Real}
    der_order::Integer
    acc_order::Integer
    num_coeffs::Integer
    num_points::Integer
    wts::Matrix{TR}
    one_over_dxn::TR
end

# TODO: benchmark this and implement a more efficient version
@inline function (op::FDMBoundOp{TR})(
    v::StridedVector{TRC}, point_idx::Integer
) where {TR<:AbstractFloat,TRC<:Union{TR,Complex{TR}}}
    return dot(@view(op.wts[:, point_idx]), v) * op.one_over_dxn
end

"""
    fdm_boundop(der_order::Integer, acc_order::Integer, dx::TR) where {TR<:AbstractFloat}

Create a shifted finite difference operator with specified derivative order and accuracy for boundary.

# Arguments
- `der_order::Integer`: The order of the derivative to approximate
- `acc_order::Integer`: The desired order of accuracy (must be even)
- `dx::TR`: Grid spacing
```
"""
function fdm_boundop(
    der_order::Integer, acc_order::Integer, dx::TR
) where {TR<:AbstractFloat}
    wts_left, wts_right = fdm_boundwts(TR, der_order, acc_order)
    num_coeffs = size(wts_left, 1)
    num_points = size(wts_left, 2)
    op_left = FDMBoundOp{TR}(
        der_order, acc_order, num_coeffs, num_points, wts_left, one(TR) / dx^der_order
    )
    op_right = FDMBoundOp{TR}(
        der_order, acc_order, num_coeffs, num_points, wts_right, one(TR) / dx^der_order
    )
    return op_left, op_right
end

struct FDMHermiteOp{TR<:Real}
    der_order::Integer
    acc_order::Integer
    num_coeffs::Integer
    num_side::Integer
    Dwts::Vector{TR}
    Ewts::Vector{TR}
    one_over_dxn::TR
    one_over_dxnm1::TR
end

# TODO: benchmark this and implement a more efficient version
function (op::FDMHermiteOp{TR})(
    f::StridedVector{TRC}, df::StridedVector{TRC}
) where {TR<:AbstractFloat,TRC<:Union{TR,Complex{TR}}}
    return dot(op.Dwts, f) * op.one_over_dxn + dot(op.Ewts, df) * op.one_over_dxnm1
end

"""
    fdm_hermiteop(der_order::Integer, acc_order::Integer, dx::TR) where {TR<:AbstractFloat}

Create a Hermite finite difference operator with specified derivative order and accuracy.

# Arguments
- `der_order::Integer`: The order of the derivative to approximate
- `acc_order::Integer`: The desired order of accuracy (must be even)
- `dx::TR`: Grid spacing
"""
function fdm_hermiteop(
    der_order::Integer, acc_order::Integer, dx::TR
) where {TR<:AbstractFloat}
    Dwts, Ewts = fdm_hermitewts(TR, der_order, acc_order)
    num_coeffs = length(Dwts)
    num_side = div(num_coeffs - 1, 2)
    return FDMHermiteOp{TR}(
        der_order,
        acc_order,
        num_coeffs,
        num_side,
        Dwts,
        Ewts,
        one(TR) / dx^der_order,
        one(TR) / dx^(der_order - 1),
    )
end

struct FDMDissOp{TR<:Real}
    diss_order::Integer
    num_coeffs::Integer
    num_side::Integer
    wts::Vector{TR}
    σ_over_dx::TR
end

# TODO: benchmark this and implement a more efficient version
@inline function (op::FDMDissOp{TR})(
    v::StridedVector{TRC}
) where {TR<:AbstractFloat,TRC<:Union{TR,Complex{TR}}}
    return dot(op.wts, v) * op.σ_over_dx
end

"""
    fdm_dissop(diss_order::Integer, σ::TR, dx::TR) where {TR<:AbstractFloat}

Create a finite difference dissipation operator with specified order.

# Arguments
- `diss_order::Integer`: The order of the dissipation operator
- `σ::TR`: Dissipation strength
- `dx::TR`: Grid spacing
"""
function fdm_dissop(diss_order::Integer, σ::TR, dx::TR) where {TR<:AbstractFloat}
    wts = fdm_disswts(diss_order)
    num_coeffs = length(wts)
    num_side = div(num_coeffs - 1, 2)
    return FDMDissOp{TR}(diss_order, num_coeffs, num_side, wts, σ / dx)
end

export FDMCentralOp, fdm_centralop
export FDMBoundOp, fdm_boundop
export FDMHermiteOp, fdm_hermiteop
export FDMDissOp, fdm_dissop
