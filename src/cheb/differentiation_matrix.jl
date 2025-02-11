@doc raw"""
    cheb_differentiation_matrix(grid::ChebyshevGrid{TF,TNode}) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode}

Compute angles for Chebyshev grid of given type.
"""
function cheb_differentiation_matrix(
    grid::ChebyshevGrid{TF,TNode}
) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode}
    return cheb_differentiation_matrix(grid.node, TF, length(grid))
end

export cheb_differentiation_matrix
