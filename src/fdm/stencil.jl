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
    return fdm_fornbergwts(der_order, zero(T), local_grid; dfdx=true)
end

function fdm_hermite(der_order::TI, acc_order::TI) where {TI<:Integer}
    return fdm_hermite(Rational{TI}, der_order, acc_order)
end

export fdm_central, fdm_hermite
