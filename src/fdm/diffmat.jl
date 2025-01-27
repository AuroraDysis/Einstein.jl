"""
    fdm_diffmat(::Type{TR}, der_order::Integer, acc_order::Integer, x_min::TR, x_max::TR, dx::TR) where {TR<:Real}

Construct a finite difference matrix for numerical differentiation.

# Arguments
- `TR`: Type parameter for the real number type to be used
- `der_order`: Order of the derivative to approximate
- `acc_order`: Order of accuracy for the approximation
- `x_min`: Minimum value of the domain
- `x_max`: Maximum value of the domain
- `dx`: Grid spacing

# Returns
- A banded matrix representing the finite difference operator with the specified derivative
  and accuracy orders. The matrix can be applied to a vector of function values to compute
  numerical derivatives.

# Notes
- Uses central difference schemes in the interior points
- Employs one-sided differences near boundaries
- The resulting matrix has dimensions n×n where n = round(Int, (x_max - x_min) / dx) + 1
"""
function fdm_diffmat(
    ::Type{TR}, der_order::Integer, acc_order::Integer, x_min::TR, x_max::TR, dx::TR;
) where {TR<:Real}
    n = round(Int, (x_max - x_min) / dx) + 1

    x_grid_end = x_min + (n - 1) * dx
    @argcheck (x_max - x_grid_end) < 10 * eps(TR) "Grid endpoint mismatch: |x_max - x_grid_end| = $(abs(x_max - x_grid_end)) exceeds tolerance ($(10 * eps(TR))). Consider adjusting dx to ensure x_max is reached precisely."

    num_coeffs = fdm_centralnum(der_order, acc_order)
    num_side = div(num_coeffs - 1, 2)
    num_boundcoeffs = fdm_boundnum(der_order, acc_order)

    diffmat = BandedMatrix(Zeros{TR}(n, n), (num_boundcoeffs, num_boundcoeffs))

    D = fdm_centralwts(TR, der_order, acc_order)
    D_left, D_right = fdm_boundwts(TR, der_order, acc_order)

    @inbounds for i in 1:num_side
        diffmat[i, 1:num_boundcoeffs] .= @view(D_left[:, i])
        diffmat[end - num_side + i, (end - num_boundcoeffs + 1):end] = D_right[:, i]
    end

    @inbounds for i in (num_side + 1):(n - num_side)
        diffmat[i, (i - num_side):(i + num_side)] .= D
    end

    return diffmat
end

export fdm_diffmat
