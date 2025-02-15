"""
    cheb_synthesis_matrix(grid::ChebyshevGrid)

Construct the synthesis matrix S that transforms Chebyshev coefficients to function values.
"""
function cheb_synthesis_matrix(
    grid::ChebyshevGrid{TF,TNode}
) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode}
    return _cheb_synthesis_matrix(grid.node, TF, length(grid))
end

export cheb_synthesis_matrix
