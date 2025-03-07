import Base: *

abstract type AbstractHermiteFiniteDifferenceOperator{TR<:Real} end
abstract type AbstractFiniteDifferenceOperator{TR<:Real} end

struct HermiteFiniteDifferenceOperator{TR<:Real,Width,HalfWidth,BoundaryWidth} <:
       AbstractHermiteFiniteDifferenceOperator{TR}
    D_weights::SVector{Width,TR}
    E_weights::SVector{Width,TR}
    D_left_weights::SMatrix{HalfWidth,BoundaryWidth,TR}
    D_right_weights::SMatrix{HalfWidth,BoundaryWidth,TR}
    E_left_weights::SMatrix{HalfWidth,BoundaryWidth,TR}
    E_right_weights::SMatrix{HalfWidth,BoundaryWidth,TR}
    D_factor::Base.RefValue{TR}
    E_factor::Base.RefValue{TR}
    derivative_order::Integer
    accuracy_order::Integer
end

struct FiniteDifferenceDerivativeOperator{TR<:Real,Width,HalfWidth,BoundaryWidth} <:
       AbstractFiniteDifferenceOperator{TR}
    weights::SVector{Width,TR}
    left_weights::SMatrix{HalfWidth,BoundaryWidth,TR}
    right_weights::SMatrix{HalfWidth,BoundaryWidth,TR}
    factor::Base.RefValue{TR}
    derivative_order::Integer
    accuracy_order::Integer
end

struct FiniteDifferenceDissipationOperator{TR<:Real,Width,HalfWidth,BoundaryWidth} <:
       AbstractFiniteDifferenceOperator{TR}
    weights::SVector{Width,TR}
    left_weights::SMatrix{HalfWidth,BoundaryWidth,TR}
    right_weights::SMatrix{HalfWidth,BoundaryWidth,TR}
    factor::Base.RefValue{TR}
    dissipation_order::Integer
end

abstract type ConvolveMode end
struct ConvolveAssign <: ConvolveMode end
struct ConvolveAdd <: ConvolveMode end

function fdm_apply_operator!(
    df::StridedArray{TR},
    op::AbstractFiniteDifferenceOperator{TR},
    f::StridedArray{TR},
    mode::Mode=ConvolveAssign(),
) where {TR<:Real,Mode<:ConvolveMode}
    (; weights, left_weights, right_weights, factor) = op
    fdm_convolve_interior!(df, f, weights, factor[], mode)
    fdm_convolve_boundary!(df, f, left_weights, right_weights, factor[], mode)
    return nothing
end

function fdm_apply_operator!(
    df::StridedArray{TR},
    op::AbstractFiniteDifferenceOperator{TR},
    f::StridedArray{TR},
    jacobian::StridedArray{TR},
    mode::Mode=ConvolveAssign(),
) where {TR<:Real,Mode<:ConvolveMode}
    (; weights, left_weights, right_weights, factor) = op
    fdm_convolve_interior!(df, f, weights, factor[], jacobian, mode)
    fdm_convolve_boundary!(df, f, left_weights, right_weights, factor[], jacobian, mode)
    return nothing
end

function fdm_apply_operator!(
    ddf::StridedArray{TR},
    op::AbstractHermiteFiniteDifferenceOperator{TR},
    f::StridedArray{TR},
    df::StridedArray{TR},
    mode::Mode=ConvolveAssign(),
) where {TR<:Real,Mode<:ConvolveMode}
    (;
        D_weights,
        E_weights,
        D_left_weights,
        E_left_weights,
        D_right_weights,
        E_right_weights,
        D_factor,
        E_factor,
    ) = op
    fdm_convolve_interior!(ddf, f, D_weights, D_factor[], mode)
    fdm_convolve_interior!(ddf, df, E_weights, E_factor[], ConvolveAdd())
    fdm_convolve_boundary!(ddf, f, D_left_weights, D_right_weights, D_factor[], mode)
    fdm_convolve_boundary!(
        ddf, df, E_left_weights, E_right_weights, E_factor[], ConvolveAdd()
    )
    return nothing
end

function fdm_apply_operator!(
    ddf::StridedArray{TR},
    op::AbstractHermiteFiniteDifferenceOperator{TR},
    f::StridedArray{TR},
    df::StridedArray{TR},
    jac_D::StridedArray{TR},
    jac_E::StridedArray{TR},
    mode::Mode=ConvolveAssign(),
) where {TR<:Real,Mode<:ConvolveMode}
    (;
        D_weights,
        E_weights,
        D_left_weights,
        E_left_weights,
        D_right_weights,
        E_right_weights,
        D_factor,
        E_factor,
    ) = op
    fdm_convolve_interior!(ddf, f, D_weights, D_factor[], jac_D, mode)
    fdm_convolve_interior!(ddf, df, E_weights, E_factor[], jac_E, ConvolveAdd())
    fdm_convolve_boundary!(
        ddf, f, D_left_weights, D_right_weights, D_factor[], jac_D, mode
    )
    fdm_convolve_boundary!(
        ddf, df, E_left_weights, E_right_weights, E_factor[], jac_E, ConvolveAdd()
    )
    return nothing
end

function *(
    op::AbstractFiniteDifferenceOperator{TR}, f::StridedArray{TR,N}
) where {TR<:Real,N}
    df = Array{TR}(undef, size(f))::Array{TR,N}
    fdm_apply_operator!(df, op, f)
    return df
end

function fdm_convolve_boundary!(
    out::StridedArray{TR},
    in::StridedArray{TR},
    left_weights::SMatrix{HalfWidth,BoundaryWidth,TR},
    right_weights::SMatrix{HalfWidth,BoundaryWidth,TR},
    factor::TR,
    ::Mode=ConvolveAssign(),
) where {TR<:Real,HalfWidth,BoundaryWidth,Mode<:ConvolveMode}
    out_left, out_right = @views out[1:HalfWidth, :], out[(end - HalfWidth + 1):end, :]
    in_left, in_right = @views in[1:BoundaryWidth, :], in[(end - BoundaryWidth + 1):end, :]
    if Mode == ConvolveAssign
        out_left .= factor .* (left_weights * in_left)
        out_right .= factor .* (right_weights * in_right)
    else # Add Mode
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
    factor::TR,
    jacobian::StridedArray{TR},
    ::Mode=ConvolveAssign(),
) where {TR<:Real,HalfWidth,BoundaryWidth,Mode<:ConvolveMode}
    out_left, out_right = @views out[1:HalfWidth, :], out[(end - HalfWidth + 1):end, :]
    in_left, in_right = @views in[1:BoundaryWidth, :], in[(end - BoundaryWidth + 1):end, :]
    jac_left = @view jacobian[1:HalfWidth, :]
    jac_right = @view jacobian[(end - HalfWidth + 1):end, :]
    if Mode == ConvolveAssign
        out_left .= factor .* jac_left .* (left_weights * in_left)
        out_right .= factor .* jac_right .* (right_weights * in_right)
    else # Add Mode
        out_left .+= factor .* jac_left .* (left_weights * in_left)
        out_right .+= factor .* jac_right .* (right_weights * in_right)
    end
    return nothing
end

@generated function fdm_convolve_interior!(
    out::StridedArray{TR},
    in::StridedArray{TR},
    weights::SVector{width,TR},
    factor::TR,
    ::Mode=ConvolveAssign(),
) where {width,TR<:Real,Mode<:ConvolveMode}
    half_width = width ÷ 2
    ex = :(weights[1] * in[begin:(end - $width + 1), :])

    for i in 2:width
        ex = :($ex + weights[$i] * in[(begin - 1 + $i):(end - $width + $i), :])
    end

    quote
        if Mode == ConvolveAssign
            @.. out[(begin + $half_width):(end - $half_width), :] = $ex * factor
        else # Add Mode
            @.. out[(begin + $half_width):(end - $half_width), :] += $ex * factor
        end
        return nothing
    end
end

@generated function fdm_convolve_interior!(
    out::StridedArray{TR},
    in::StridedArray{TR},
    weights::SVector{width,TR},
    factor::TR,
    jacobian::StridedArray{TR},
    ::Mode=ConvolveAssign(),
) where {width,TR<:Real,Mode<:ConvolveMode}
    half_width = width ÷ 2
    ex = :(weights[1] * in[begin:(end - $width + 1), :])

    for i in 2:width
        ex = :($ex + weights[$i] * in[(begin - 1 + $i):(end - $width + $i), :])
    end

    quote
        if Mode == ConvolveAssign
            @.. broadcast = true out[(begin + $half_width):(end - $half_width), :] =
                factor * jacobian[(begin + $half_width):(end - $half_width)] * $ex
        else # Add Mode
            @.. broadcast = true out[(begin + $half_width):(end - $half_width), :] +=
                factor * jacobian[(begin + $half_width):(end - $half_width)] * $ex
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

function fdm_hermite_derivative_operator(
    ::Type{TR}, derivative_order::Integer, accuracy_order::Integer, dx::TR
) where {TR<:Real}
    D_weights, E_weights = fdm_hermite_weights(TR, derivative_order, accuracy_order)
    D_left_weights, E_left_weights, D_right_weights, E_right_weights = fdm_hermite_boundary_weights(
        TR, derivative_order, accuracy_order
    )

    width = length(D_weights)
    half_width = div(width - 1, 2)
    boundary_width = size(D_left_weights, 2)
    D_factor = inv(dx^derivative_order)
    E_factor = inv(dx^(derivative_order - 1))
    return HermiteFiniteDifferenceOperator{TR,width,half_width,boundary_width}(
        D_weights,
        E_weights,
        D_left_weights,
        D_right_weights,
        E_left_weights,
        E_right_weights,
        Ref(D_factor),
        Ref(E_factor),
        derivative_order,
        accuracy_order,
    )
end

function fdm_hermite_derivative_operator(
    derivative_order::Integer, accuracy_order::Integer, dx::Float64
)
    return fdm_hermite_derivative_operator(Float64, derivative_order, accuracy_order, dx)
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
    fdm_operator_matrix(op::AbstractFiniteDifferenceOperator{TR}; boundary::Bool=false, transpose::Bool=false) -> BandedMatrix{TR}

Create a banded matrix representation of the finite difference operator.

# Arguments
- `op::AbstractFiniteDifferenceOperator{TR}`: The finite difference operator
- `boundary::Bool=true`: Whether to include boundary weights
- `transpose::Bool=false`: Whether to transpose the matrix
"""
function fdm_operator_matrix(
    op::AbstractFiniteDifferenceOperator{TR},
    n::Integer,
    boundary::Bool=true,
    transpose::Bool=false,
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
    HermiteFiniteDifferenceOperator,
    FiniteDifferenceDissipationOperator,
    fdm_derivative_operator,
    fdm_hermite_derivative_operator,
    fdm_dissipation_operator,
    fdm_operator_matrix,
    ConvolveAdd,
    ConvolveAssign,
    fdm_apply_operator!
export fdm_convolve_boundary!, fdm_convolve_interior!
