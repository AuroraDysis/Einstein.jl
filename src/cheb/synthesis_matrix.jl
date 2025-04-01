"""
    cheb_synthesis_matrix(grid::ChebyshevGrid)

Construct the synthesis matrix S that transforms Chebyshev coefficients to function values.
"""
function cheb_synthesis_matrix(
    grid::ChebyshevGrid{TF,Basis}
) where {TF<:AbstractFloat,Basis<:AbstractBasisFunction}
    return _cheb_synthesis_matrix(grid.basis, TF, length(grid))
end

export cheb_synthesis_matrix
