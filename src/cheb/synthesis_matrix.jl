function cheb_synthesis_matrix(
    grid::ChebyshevGrid{TF,TNode}
) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode}
    return _cheb_synthesis_matrix(grid.node, TF, length(grid))
end
