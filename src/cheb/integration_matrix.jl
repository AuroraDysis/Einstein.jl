@doc raw"""
    cheb_integration_matrix(grid::ChebyshevGrid{TF,TNode}) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode}

Compute the Chebyshev integration matrix that maps function values at `n` Chebyshev points to values of the integral of the interpolating polynomial at those points, with the convention that the first value is zero.
"""
function cheb_integration_matrix(
    grid::ChebyshevGrid{TF,TNode}
) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode}
    return _cheb_integration_matrix(
        grid.node, TF, length(grid), grid.lower_bound, grid.upper_bound
    )
end

export cheb_integration_matrix
