"""
    fdm_central_width(derivative_order::Integer, accuracy_order::Integer)

Calculate the number of coefficients needed for central FDM stencil.

# Arguments
- `derivative_order::Integer`: Order of the derivative
- `accuracy_order::Integer`: Order of accuracy

# References
- [Finite difference coefficient - Wikipedia](https://en.wikipedia.org/wiki/Finite_difference_coefficient)
"""
function fdm_central_width(derivative_order::Integer, accuracy_order::Integer)
    @boundscheck begin
        @argcheck derivative_order >= 1 "Derivative order must be at least 1"
        @argcheck accuracy_order % 2 == 0 "Only even orders are supported for central FDM stencils."
    end

    # https://github.com/maroba/findiff/blob/master/findiff/coefs.py
    # coefficients number for central difference
    return 2 * div(derivative_order + 1, 2) - 1 + accuracy_order
end

"""
    fdm_hermite_width(derivative_order::Integer, accuracy_order::Integer)

Calculate the number of coefficients needed for Hermite FDM stencil.

# Arguments
- `derivative_order::Integer`: Order of the derivative
- `accuracy_order::Integer`: Order of accuracy

# References
- [fornberg2021algorithm](@citet*)
"""
function fdm_hermite_width(derivative_order::Integer, accuracy_order::Integer)
    @boundscheck begin
        @argcheck derivative_order >= 2 "Only derivative order greater than or equal to 2 are supported for Hermite-type finite difference."

        if mod(div(derivative_order, 2), 2) == 1
            # accuracy_order must be 4,8,12... for der order 2,3,6,7,10,11...
            @argcheck accuracy_order % 4 == 0 "Only accuracy_order % 4 == 0 are supported for Hermite-type finite difference with der order 2,3,6,7,10,11..."
        else
            # accuracy_order must be 2,6,10... for der order 4,5,8,9,12...
            @argcheck accuracy_order % 4 == 2 "Only accuracy_order % 4 == 2 are supported for Hermite-type finite difference with der order 4,5,8,9,12..."
        end
    end

    return div(derivative_order, 2) + div(accuracy_order, 2)
end

"""
    fdm_boundary_width(derivative_order::Integer, accuracy_order::Integer)

Calculate the number of coefficients needed for shifted boundary FDM stencil.

# Arguments
- `derivative_order::Integer`: Order of the derivative
- `accuracy_order::Integer`: Order of accuracy
"""
function fdm_boundary_width(derivative_order::Integer, accuracy_order::Integer)
    @boundscheck begin
        @argcheck derivative_order >= 1 "Derivative order must be at least 1"
        @argcheck accuracy_order >= 1 "Accuracy order must be at least 1"
    end

    return derivative_order + accuracy_order
end

"""
    fdm_dissipation_width(dissipation_order::Integer)

Calculate the number of coefficients needed for Kreiss-Oliger dissipation (interior).
"""
function fdm_dissipation_width(dissipation_order::Integer)
    return fdm_central_width(dissipation_order, 2)
end

"""
    fdm_dissipation_boundary_width(dissipation_order::Integer)

Calculate the number of coefficients needed for Kreiss-Oliger dissipation (boundary).
"""
function fdm_dissipation_boundary_width(dissipation_order::Integer)
    return fdm_boundary_width(dissipation_order, 2)
end

export fdm_central_width, fdm_hermite_width, fdm_boundary_width, fdm_dissipation_width, fdm_dissipation_boundary_width
