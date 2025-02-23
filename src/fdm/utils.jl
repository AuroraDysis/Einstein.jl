"""
    fdm_central_width(der_order::Integer, acc_order::Integer)

Calculate the number of coefficients needed for central FDM stencil.

# Arguments
- `der_order::Integer`: Order of the derivative
- `acc_order::Integer`: Order of accuracy

# References
- [Finite difference coefficient - Wikipedia](https://en.wikipedia.org/wiki/Finite_difference_coefficient)
"""
function fdm_central_width(der_order::Integer, acc_order::Integer)
    # https://github.com/maroba/findiff/blob/master/findiff/coefs.py
    # coefficients number for central difference
    return 2 * div(der_order + 1, 2) - 1 + acc_order
end

"""
    fdm_hermite_width(der_order::Integer, acc_order::Integer)

Calculate the number of coefficients needed for Hermite FDM stencil.

# Arguments
- `der_order::Integer`: Order of the derivative
- `acc_order::Integer`: Order of accuracy

# References
- [fornberg2021algorithm](@citet*)
"""
function fdm_hermite_width(der_order::Integer, acc_order::Integer)
    return return div(der_order, 2) + div(acc_order, 2)
end

"""
    fdm_boundary_width(der_order::Integer, acc_order::Integer)

Calculate the number of coefficients needed for shifted boundary FDM stencil.

# Arguments
- `der_order::Integer`: Order of the derivative
- `acc_order::Integer`: Order of accuracy
"""
function fdm_boundary_width(der_order::Integer, acc_order::Integer)
    return der_order + acc_order
end

"""
    fdm_dissipation_width(diss_order::Integer)

Calculate the number of coefficients needed for Kreiss-Oliger dissipation (interior).
"""
function fdm_dissipation_width(diss_order::Integer)
    return fdm_central_width(diss_order, 2)
end

"""
    fdm_dissipation_boundary_width(diss_order::Integer)

Calculate the number of coefficients needed for Kreiss-Oliger dissipation (boundary).
"""
function fdm_dissipation_boundary_width(diss_order::Integer)
    return fdm_boundary_width(diss_order, 2)
end

export fdm_central_width, fdm_hermite_width, fdm_boundary_width, fdm_dissipation_width, fdm_dissipation_boundary_width
