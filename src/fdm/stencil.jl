"""
    fdm_central([T=Rational{TI}], der_order::TI, acc_order::TI) where {T<:Real, TI<:Integer}

Generate central finite difference coefficients for a given derivative and accuracy order.

# Arguments
- `der_order::Integer`: The order of the derivative to approximate
- `acc_order::Integer`: The desired order of accuracy (must be even)

# Returns
Vector of rational coefficients for the finite difference stencil
"""
function fdm_central(::Type{T}, der_order::Integer, acc_order::Integer) where {T<:Real}
    @argcheck acc_order % 2 == 0 "Only even orders are supported for central FDM stencils."

    num_coeffs = fdm_centralnum(der_order, acc_order)
    num_side = div(num_coeffs - 1, 2)
    local_grid = collect(T, (-num_side):num_side)
    return fdm_fornbergwts(der_order, zero(T), local_grid)
end

function fdm_central(der_order::TI, acc_order::TI) where {TI<:Integer}
    return fdm_central(Rational{TI}, der_order, acc_order)
end

"""
    fdm_hermite([T=Rational{TI}], der_order::TI, acc_order::TI) where {T<:Real, TI<:Integer}

Generate Hermite-type finite difference coefficients that include function value and derivative information.

# Arguments
- `der_order::Integer`: The order of the derivative to approximate (must be ≥ 2)
- `acc_order::Integer`: The desired order of accuracy
    * For der_order 2,3,6,7,10,11...: acc_order must be 4,8,12...
    * For der_order 4,5,8,9,12...: acc_order must be 2,6,10...

# Returns
Vector of rational coefficients for the Hermite-type finite difference stencil
"""
function fdm_hermite(::Type{T}, der_order::Integer, acc_order::Integer) where {T<:Real}
    @argcheck der_order >= 2 "Only derivative order greater than or equal to 2 are supported for Hermite-type finite difference."

    if mod(div(der_order, 2), 2) == 1
        # acc_order must be 4,8,12... for der order 2,3,6,7,10,11...
        @argcheck acc_order % 4 == 0 "Only acc_order % 4 == 0 are supported for Hermite-type finite difference with der order 2,3,6,7,10,11..."
    else
        # acc_order must be 2,6,10... for der order 4,5,8,9,12...
        @argcheck acc_order % 4 == 2 "Only acc_order % 4 == 2 are supported for Hermite-type finite difference with der order 4,5,8,9,12..."
    end

    num_coeffs = fdm_hermitenum(der_order, acc_order)
    num_side = div(num_coeffs - 1, 2)
    local_grid = collect(T, (-num_side):num_side)
    return fdm_fornbergwts(der_order, zero(T), local_grid; hermite=true)
end

function fdm_hermite(der_order::TI, acc_order::TI) where {TI<:Integer}
    return fdm_hermite(Rational{TI}, der_order, acc_order)
end

"""
    fdm_extrapwts_left(extrap_order::Int)

Generate weights for left-sided extrapolation of order `extrap_order`.

# Arguments
- `extrap_order::Int`: Order of extrapolation

# Returns
Vector of rational coefficients for left-sided extrapolation
"""
function fdm_extrapwts_left(extrap_order::Int)
    return fdm_fornbergwts(0, 0, 1:extrap_order)
end

"""
    fdm_extrapwts_right(extrap_order::Int)

Generate weights for right-sided extrapolation of order `extrap_order`.

# Arguments
- `extrap_order::Int`: Order of extrapolation

# Returns
Vector of rational coefficients for right-sided extrapolation
"""
function fdm_extrapwts_right(extrap_order::Int)
    return fdm_fornbergwts(0, 0, extrap_order:-1:1)
end

"""
    fdm_boundwts([T=Rational{TI}], der_order::TI, acc_order::TI) where {T<:Real, TI<:Integer}

Generate finite difference coefficients for shifted boundary conditions.

# Arguments
- `der_order::Integer`: The order of the derivative to approximate
- `acc_order::Integer`: The desired order of accuracy

# Returns
Tuple of left and right shifted boundary finite difference coefficients
The coefficients are stored in a matrix with the columns representing the different grid points.
The columns are ordered from the leftmost grid point to the rightmost grid point.
"""
function fdm_boundwts(::Type{T}, der_order::Integer, acc_order::Integer) where {T<:Real}
    num_coeffs = fdm_boundnum(der_order, acc_order)
    num_central = fdm_centralnum(der_order, acc_order)
    num_side = div(num_central - 1, 2)

    D_left = zeros(T, num_coeffs, num_side)
    D_right = zeros(T, num_coeffs, num_side)

    local_grid = zeros(T, num_coeffs)
    @inbounds for i in 1:num_side
        for j in 1:num_coeffs
            local_grid[j] = -(i - 1) + (j - 1)
        end
        D_left[:, i] = fdm_fornbergwts(der_order, zero(T), local_grid)

        for j in 1:num_coeffs
            local_grid[end - j + 1] = (i - 1) - (j - 1)
        end
        D_right[:, end - i + 1] = fdm_fornbergwts(der_order, zero(T), local_grid)
    end

    return D_left, D_right
end

function fdm_boundwts(der_order::TI, acc_order::TI) where {TI<:Integer}
    return fdm_boundwts(Rational{TI}, der_order, acc_order)
end

export fdm_central, fdm_hermite, fdm_extrapwts_right, fdm_extrapwts_left, fdm_boundwts
