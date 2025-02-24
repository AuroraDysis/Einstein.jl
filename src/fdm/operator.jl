import Base: *
using Static

abstract type FiniteDifferenceOperator{TR<:Real} end

struct FiniteDifferenceDerivativeOperator{TR<:Real,Width,HalfWidth,BoundaryWidth} <:
       FiniteDifferenceOperator{TR}
    weights::SVector{Width,TR}
    left_weights::SMatrix{HalfWidth,BoundaryWidth,TR}
    right_weights::SMatrix{HalfWidth,BoundaryWidth,TR}
    factor::Base.RefValue{TR}
    derivative_order::Integer
    accuracy_order::Integer
end

struct FiniteDifferenceDissipationOperator{TR<:Real,Width,HalfWidth,BoundaryWidth} <:
       FiniteDifferenceOperator{TR}
    weights::SVector{Width,TR}
    left_weights::SMatrix{HalfWidth,BoundaryWidth,TR}
    right_weights::SMatrix{HalfWidth,BoundaryWidth,TR}
    factor::Base.RefValue{TR}
    dissipation_order::Integer
end

const Add = StaticInt{1}
const Assign = StaticInt{0}

function mul!(
    df::StridedArray{TR}, op::FiniteDifferenceOperator{TR}, f::StridedArray{TR}
) where {TR<:Real}
    (; weights, left_weights, right_weights, factor) = op
    fdm_convolve_interior!(df, f, weights, factor[])
    fdm_convolve_boundary!(df, f, left_weights, right_weights, factor[])
    return nothing
end

function *(op::FiniteDifferenceOperator{TR}, f::StridedArray{TR,N}) where {TR<:Real,N}
    df = Array{TR}(undef, size(f))::Array{TR,N}
    mul!(df, op, f)
    return df
end

function fdm_convolve_boundary!(
    out::StridedArray{TR},
    in::StridedArray{TR},
    left_weights::SMatrix{HalfWidth,BoundaryWidth,TR},
    right_weights::SMatrix{HalfWidth,BoundaryWidth,TR},
    factor::TR,
    mode::StaticInt=Assign(),
) where {TR<:Real,HalfWidth,BoundaryWidth}
    out_left, out_right = @views out[1:HalfWidth, :], out[(end - HalfWidth + 1):end, :]
    in_left, in_right = @views in[1:BoundaryWidth, :], in[(end - BoundaryWidth + 1):end, :]
    if mode == Assign()
        out_left .= factor .* (left_weights * in_left)
        out_right .= factor .* (right_weights * in_right)
    else # Add mode
        out_left .+= factor .* (left_weights * in_left)
        out_right .+= factor .* (right_weights * in_right)
    end
    return nothing
end

function fdm_convolve_boundary!(
    out::StridedArray{TR},
    in::StridedArray{TR},
    left_weights::SMatrix{HalfWidth,BoundaryWidth,TR},
    right_weights::SMatrix{HalfWidth,BoundaryWidth,TR},
    factors::StridedVector{TR},
    mode::StaticInt=Assign(),
) where {TR<:Real,HalfWidth,BoundaryWidth}
    out_left, out_right = @views out[1:HalfWidth, :], out[(end - HalfWidth + 1):end, :]
    in_left, in_right = @views in[1:BoundaryWidth, :], in[(end - BoundaryWidth + 1):end, :]
    factors_left, factors_right = @views factors[1:HalfWidth],
    factors[(end - HalfWidth + 1):end]
    if mode == Assign()
        out_left .= factors_left .* (left_weights * in_left)
        out_right .= factors_right .* (right_weights * in_right)
    else # Add mode
        out_left .+= factors_left .* (left_weights * in_left)
        out_right .+= factors_right .* (right_weights * in_right)
    end
    return nothing
end

@generated function fdm_convolve_interior!(
    out::StridedArray{TR},
    in::StridedArray{TR},
    weights::SVector{width,TR},
    factor::TR,
    mode::StaticInt=Assign(),
) where {width,TR<:Real}
    half_width = width ÷ 2
    ex = :(weights[1] * in[begin:(end - $width + 1), :])

    for i in 2:width
        ex = :($ex + weights[$i] * in[(begin - 1 + $i):(end - $width + $i), :])
    end

    quote
        if mode == Assign()
            @.. out[(begin + $half_width):(end - $half_width), :] = $ex * factor
        else # Add mode
            @.. out[(begin + $half_width):(end - $half_width), :] += $ex * factor
        end
        return nothing
    end
end

@generated function fdm_convolve_interior!(
    out::StridedArray{TR},
    in::StridedArray{TR},
    weights::SVector{width,TR},
    factors::StridedVector{TR},
    mode::StaticInt=Assign(),
) where {width,TR<:Real}
    half_width = width ÷ 2
    ex = :(weights[1] * in[begin:(end - $width + 1), :])

    for i in 2:width
        ex = :($ex + weights[$i] * in[(begin - 1 + $i):(end - $width + $i), :])
    end

    quote
        if mode == Assign()
            @.. broadcast = true out[(begin + $half_width):(end - $half_width), :] =
                factors[(begin + $half_width):(end - $half_width)] * $ex
        else # Add mode
            @.. broadcast = true out[(begin + $half_width):(end - $half_width), :] +=
                factors[(begin + $half_width):(end - $half_width)] * $ex
        end
        return nothing
    end
end

"""
    fdm_derivative_operator([TR=Float64], derivative_order::Integer, accuracy_order::Integer, dx::TR) -> FiniteDifferenceDerivativeOperator{TR}

Create a finite difference derivative operator with specified derivative and accuracy orders.

# Arguments
- `TR`: The element type of the operator
- `derivative_order::Integer`: The order of the derivative
- `accuracy_order::Integer`: The order of accuracy
- `dx::TR`: The grid spacing
"""
function fdm_derivative_operator(
    ::Type{TR}, derivative_order::Integer, accuracy_order::Integer, dx::TR
) where {TR<:Real}
    weights = fdm_central_weights(TR, derivative_order, accuracy_order)
    left_weights, right_weights = fdm_boundary_weights(TR, derivative_order, accuracy_order)
    width = length(weights)
    half_width = div(width - 1, 2)
    boundary_width = size(left_weights, 2)
    factor = inv(dx^derivative_order)
    return FiniteDifferenceDerivativeOperator{TR,width,half_width,boundary_width}(
        SVector{width,TR}(weights),
        SMatrix{half_width,boundary_width,TR}(left_weights),
        SMatrix{half_width,boundary_width,TR}(right_weights),
        Ref(factor),
        derivative_order,
        accuracy_order,
    )
end

function fdm_derivative_operator(
    derivative_order::Integer, accuracy_order::Integer, dx::Float64
)
    return fdm_derivative_operator(Float64, derivative_order, accuracy_order, dx)
end

"""
    fdm_dissipation_operator([TR=Float64], dissipation_order::Integer, σ::TR, dx::TR) -> FiniteDifferenceDerivativeOperator{TR}

Create a finite difference dissipation operator with specified dissipation order.

# Arguments
- `TR`: The element type of the operator
- `dissipation_order::Integer`: The order of the dissipation operator
- `σ::TR`: The dissipation coefficient
- `dx::TR`: The grid spacing
"""
function fdm_dissipation_operator(
    ::Type{TR}, dissipation_order::Integer, σ::TR, dx::TR
) where {TR<:Real}
    weights = fdm_dissipation_weights(TR, dissipation_order)
    left_weights, right_weights = fdm_dissipation_boundary_weights(TR, dissipation_order)
    width = length(weights)
    half_width = div(width - 1, 2)
    boundary_width = size(left_weights, 2)
    factor = σ / dx
    return FiniteDifferenceDissipationOperator{TR,width,half_width,boundary_width}(
        SVector{width,TR}(weights),
        SMatrix{half_width,boundary_width,TR}(left_weights),
        SMatrix{half_width,boundary_width,TR}(right_weights),
        Ref(factor),
        dissipation_order,
    )
end

function fdm_dissipation_operator(dissipation_order::Integer, σ::Float64, dx::Float64)
    return fdm_dissipation_operator(Float64, dissipation_order, σ, dx)
end

"""
    fdm_operator_matrix(op::FiniteDifferenceOperator{TR}; boundary::Bool=false, transpose::Bool=false) -> BandedMatrix{TR}

Create a banded matrix representation of the finite difference operator.

# Arguments
- `op::FiniteDifferenceOperator{TR}`: The finite difference operator
- `boundary::Bool=true`: Whether to include boundary weights
- `transpose::Bool=false`: Whether to transpose the matrix
"""
function fdm_operator_matrix(
    op::FiniteDifferenceOperator{TR}, n::Integer, boundary::Bool=true, transpose::Bool=false
) where {TR<:Real}
    (; weights, left_weights, right_weights, factor) = op
    half_width = length(weights) ÷ 2
    weights *= factor[]
    left_weights *= factor[]
    right_weights *= factor[]

    if boundary
        boundary_width = size(left_weights, 2)
        mat = BandedMatrix(Zeros{TR}(n, n), (boundary_width, boundary_width))

        @inbounds for i in 1:half_width
            mat[i, 1:boundary_width] .= @view(left_weights[i, :])
            mat[end - half_width + i, (end - boundary_width + 1):end] = @view(
                right_weights[i, :]
            )
        end
    else
        mat = BandedMatrix(Zeros{TR}(n, n), (half_width, half_width))
    end

    @inbounds for i in (half_width + 1):(n - half_width)
        mat[i, (i - half_width):(i + half_width)] .= weights
    end

    if transpose
        diffmat_T = similar(mat)
        transpose!(diffmat_T, mat)
        return diffmat_T
    else
        return mat
    end
end

export FiniteDifferenceDerivativeOperator,
    FiniteDifferenceDissipationOperator,
    fdm_derivative_operator,
    fdm_dissipation_operator,
    fdm_operator_matrix,
    Add,
    Assign
export fdm_convolve_boundary!, fdm_convolve_interior!
