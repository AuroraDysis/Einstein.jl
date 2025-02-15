"""
    cheb_analysis_matrix(grid::ChebyshevGrid)

Construct the analysis matrix A that transforms function values to Chebyshev coefficients.
"""
function cheb_analysis_matrix(
    grid::ChebyshevGrid{TF,TNode}
) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode}
    return _cheb_analysis_matrix(grid.node, TF, length(grid))
end

export cheb_analysis_matrix
