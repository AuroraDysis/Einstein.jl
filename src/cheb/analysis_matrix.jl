function cheb_analysis_matrix(
    grid::ChebyshevGrid{TF,TNode}
) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode}
    return cheb_analysis_matrix(grid.node, TF, length(grid))
end
