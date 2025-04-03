"""
    fdm_dissipation_order(accuracy_order::Integer)

Calculate the order of dissipation needed for a given finite difference accuracy order [Babiuc:2007vr](@cite).
For a scheme of accuracy order 2r-2, returns dissipation order 2r.
"""
function fdm_dissipation_order(accuracy_order::Integer)
    @argcheck iseven(accuracy_order) "Only even orders are supported."
    r = div(accuracy_order + 2, 2)
    return 2r
end

"""
    fdm_dissipation_weights([TR=Rational{Int}], dissipation_order::Integer)

Calculate the weights for Kreiss-Oliger dissipation of given order [Babiuc:2007vr](@cite).
"""
function fdm_dissipation_weights(::Type{TR}, dissipation_order::Integer) where {TR<:Real}
    @argcheck iseven(dissipation_order) "Only even orders are supported."
    r = div(dissipation_order, 2)
    weights = fdm_central_weights(TR, dissipation_order, 2)
    factor = (-1)^(r + 1) / 2^(2 * r)
    @.. weights *= factor
    return weights
end

function fdm_dissipation_weights(dissipation_order::TI) where {TI<:Integer}
    return fdm_dissipation_weights(Rational{TI}, dissipation_order)
end

"""
    fdm_dissipation_boundary_weights([TR=Rational{Int}], dissipation_order::Integer)

Calculate the weights for Kreiss-Oliger dissipation of given order at the boundary [Babiuc:2007vr](@cite).
"""
function fdm_dissipation_boundary_weights(
    ::Type{TR}, dissipation_order::Integer
) where {TR<:Real}
    @argcheck iseven(dissipation_order) "Only even orders are supported."
    r = div(dissipation_order, 2)
    weights_left, weights_right = fdm_boundary_weights(TR, dissipation_order, 2)
    factor = (-1)^(r + 1) / 2^(2 * r)
    @.. weights_left *= factor
    @.. weights_right *= factor
    return weights_left, weights_right
end

function fdm_dissipation_boundary_weights(dissipation_order::TI) where {TI<:Integer}
    return fdm_dissipation_boundary_weights(Rational{TI}, dissipation_order)
end

export fdm_dissipation_order, fdm_dissipation_weights, fdm_dissipation_boundary_weights
