@doc raw"""
    cheb_angles(grid::ChebyshevGrid{TF,TNode}) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode}

Compute angles for Chebyshev grid of given type.
"""
function cheb_angles(
    grid::ChebyshevGrid{TF,TNode}
) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode}
    return _cheb_angles(grid.node, TF, length(grid))
end

export cheb_angles
