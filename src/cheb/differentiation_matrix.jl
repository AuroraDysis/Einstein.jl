@doc raw"""
    cheb_differentiation_matrix(grid::ChebyshevGrid{TF,TNode}) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode}

Compute the Chebyshev differentiation matrix that maps function values at `n` Chebyshev points to values of the derivative of the interpolating polynomial at those points.
"""
function cheb_differentiation_matrix(
    grid::ChebyshevGrid{TF,TNode}
) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode}
    return cheb_differentiation_matrix(grid.node, TF, length(grid))
end

export cheb_differentiation_matrix
