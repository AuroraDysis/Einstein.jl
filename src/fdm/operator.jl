struct FiniteDifferenceOperator{TF<:AbstractFloat,Width,HalfWidth,BoundaryWidth}
    weights::SVector{Width,TF}
    left_weights::SMatrix{HalfWidth,BoundaryWidth,TF}
    right_weights::SMatrix{HalfWidth,BoundaryWidth,TF}
    factor::Base.RefValue{TF}
end

function mul!(
    df::StridedArray{TF}, op::FiniteDifferenceOperator{TF}, f::StridedArray{TF}
) where {TF<:AbstractFloat}
    (; weights, left_weights, right_weights, factor) = op
    fdm_convolve_interior!(df, f, weights, factor[])
    fdm_convolve_boundary!(df, f, left_weights, right_weights, factor[])
    return nothing
end

function fdm_convolve_boundary!(
    out::StridedArray{TF},
    in::StridedArray{TF},
    left_weights::SMatrix{HalfWidth,BoundaryWidth,TF},
    right_weights::SMatrix{HalfWidth,BoundaryWidth,TF},
    factor::TF,
) where {TF<:AbstractFloat,HalfWidth,BoundaryWidth}
    out_left, out_right = @views out[1:HalfWidth, :], out[(end - HalfWidth + 1):end, :]
    in_left, in_right = @views in[1:BoundaryWidth, :], in[(end - BoundaryWidth + 1):end, :]
    out_left .= factor * left_weights * in_left
    out_right .= factor * right_weights * in_right
    return nothing
end

@generated function fdm_convolve_interior!(
    out::StridedArray{TF}, in::StridedArray{TF}, weights::SVector{width,TF}, factor::TF
) where {width,TF<:AbstractFloat}
    half_width = width ÷ 2
    ex = :(weights[1] * in[begin:(end - $width + 1), :])

    for i in 2:width
        ex = :($ex + weights[$i] * in[(begin - 1 + $i):(end - $width + $i), :])
    end

    quote
        @.. out[(begin + $half_width):(end - $half_width), :] = $ex * factor
        return nothing
    end
end

"""
    fdm_operator(::Type{TF}, derivative_order::Integer, accuracy_order::Integer, dx::TF) -> FiniteDifferenceOperator{TF}

Create a finite difference operator with specified derivative and accuracy orders.

# Arguments
- `TF`: The element type of the operator
- `derivative_order::Integer`: The order of the derivative
- `accuracy_order::Integer`: The order of accuracy
- `dx::TF`: The grid spacing
"""
function fdm_operator(
    ::Type{TF}, derivative_order::Integer, accuracy_order::Integer, dx::TF
) where {TF<:AbstractFloat}
    weights = fdm_central_weights(TF, derivative_order, accuracy_order)
    left_weights, right_weights = fdm_boundary_weights(TF, derivative_order, accuracy_order)
    width = length(weights)
    half_width = div(width - 1, 2)
    boundary_width = size(left_weights, 2)
    factor = inv(dx^derivative_order)
    return FiniteDifferenceOperator{TF,width,half_width,boundary_width}(
        SVector{width,TF}(weights),
        SMatrix{half_width,boundary_width,TF}(left_weights),
        SMatrix{half_width,boundary_width,TF}(right_weights),
        Ref(factor),
    )
end

"""
    fdm_dissipation_operator(::Type{TF}, dissipation_order::Integer, σ::TF, dx::TF) -> FiniteDifferenceOperator{TF}

Create a finite difference dissipation operator with specified dissipation order.

# Arguments
- `TF`: The element type of the operator
- `dissipation_order::Integer`: The order of the dissipation operator
- `σ::TF`: The dissipation coefficient
- `dx::TF`: The grid spacing
"""
function fdm_dissipation_operator(
    ::Type{TF}, dissipation_order::Integer, σ::TF, dx::TF
) where {TF<:AbstractFloat}
    weights = fdm_dissipation_weights(TF, dissipation_order)
    left_weights, right_weights = fdm_dissipation_boundary_weights(TF, dissipation_order)
    width = length(weights)
    half_width = div(width - 1, 2)
    boundary_width = size(left_weights, 2)
    factor = σ / dx
    return new{TF,width,half_width,0}(
        SVector{width,TF}(weights),
        SMatrix{half_width,boundary_width,TF}(left_weights),
        SMatrix{half_width,boundary_width,TF}(right_weights),
        Ref(factor),
    )
end

export FiniteDifferenceOperator, fdm_operator, fdm_dissipation_operator
export fdm_convolve_boundary!, fdm_convolve_interior!
