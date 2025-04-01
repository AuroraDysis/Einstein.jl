"""
    cheb_quadrature_weights(grid::ChebyshevGrid)

Compute the quadrature weights for a Chebyshev grid.
"""
function cheb_quadrature_weights(
    grid::ChebyshevGrid{TF,Basis}
) where {TF<:AbstractFloat,Basis<:AbstractBasisFunction}
    return _cheb_quadrature_weights(
        grid.basis, TF, length(grid), grid.lower_bound, grid.upper_bound
    )
end

export cheb_quadrature_weights
