"""
    fdm_diffmat(::Type{TR}, der_order::Integer, acc_order::Integer, n::Integer, with_bound::Bool=false)

Construct a finite difference matrix for numerical differentiation.

# Arguments
- `TR`: Type parameter for the real number type to be used
- `der_order`: Order of the derivative to approximate
- `acc_order`: Order of accuracy for the approximation
- `n`: Number of grid points
- `with_bound`: Flag to indicate if the matrix should include shifted boundary finite difference coefficients (default: `false`)
- `transpose`: Flag to indicate if the matrix should be transposed (default: `false`)

# Returns
- A banded matrix representing the finite difference operator with the specified derivative
  and accuracy orders. The matrix can be applied to a vector of function values to compute
  numerical derivatives.

# Notes
- Uses central difference schemes in the interior points
- Employs one-sided differences near boundaries
- The resulting matrix has dimensions n×n where n = round(Int, (upper_bound - lower_bound) / dx) + 1
"""
function fdm_diffmat(
    ::Type{TR},
    der_order::Integer,
    acc_order::Integer,
    n::Integer;
    with_bound::Bool=true,
    transpose::Bool=false,
) where {TR<:Real}
    if with_bound
        diffmat = fdm_diffmat_bound(TR, der_order, acc_order, n)
    else
        diffmat = fdm_diffmat_central(TR, der_order, acc_order, n)
    end

    if transpose
        diffmat_T = similar(diffmat)
        transpose!(diffmat_T, diffmat)
        return diffmat_T
    else
        return diffmat
    end
end

function fdm_diffmat_central(
    ::Type{TR}, der_order::Integer, acc_order::Integer, n::Integer
) where {TR<:Real}
    op = fdm_centralop(der_order, acc_order, one(TR))
    num_side = op.num_side
    wts = op.wts

    diffmat = BandedMatrix(Zeros{TR}(n, n), (num_side, num_side))
    wts = fdm_central_weights(TR, der_order, acc_order)

    @inbounds for i in (num_side + 1):(n - num_side)
        diffmat[i, (i - num_side):(i + num_side)] .= wts
    end

    return diffmat
end

function fdm_diffmat_bound(
    ::Type{TR}, der_order::Integer, acc_order::Integer, n::Integer
) where {TR<:Real}
    op = fdm_centralop(der_order, acc_order, one(TR))
    num_side = op.num_side
    wts = op.wts

    op_left, op_right = fdm_boundop(der_order, acc_order, one(TR))
    num_boundcoeffs = op_left.num_coeffs
    wts_left = op_left.wts
    wts_right = op_right.wts

    diffmat = BandedMatrix(Zeros{TR}(n, n), (num_boundcoeffs, num_boundcoeffs))

    wts = fdm_central_weights(TR, der_order, acc_order)
    wts_left, wts_right = fdm_boundary_weights(TR, der_order, acc_order)

    @inbounds for i in 1:num_side
        diffmat[i, 1:num_boundcoeffs] .= @view(wts_left[:, i])
        diffmat[end - num_side + i, (end - num_boundcoeffs + 1):end] = wts_right[:, i]
    end

    @inbounds for i in (num_side + 1):(n - num_side)
        diffmat[i, (i - num_side):(i + num_side)] .= wts
    end

    return diffmat
end

export fdm_diffmat, fdm_diffmat_central, fdm_diffmat_bound
