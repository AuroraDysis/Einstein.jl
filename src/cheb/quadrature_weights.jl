"""
    cheb_quadrature_weights(grid::ChebyshevGrid)

Compute the quadrature weights for a Chebyshev grid.
"""
function cheb_quadrature_weights(
    grid::ChebyshevGrid{TF,TNode}
) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode}
    return _cheb_quadrature_weights(
        grid.node, TF, length(grid), grid.lower_bound, grid.upper_bound
    )
end

export cheb_integration_matrix
