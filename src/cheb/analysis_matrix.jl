"""
    cheb_analysis_matrix(grid::ChebyshevGrid)

Construct the analysis matrix A that transforms function values to Chebyshev coefficients.
"""
function cheb_analysis_matrix(
    grid::ChebyshevGrid{TF,Basis}
) where {TF<:AbstractFloat,Basis<:AbstractBasisFunction}
    return _cheb_analysis_matrix(grid.basis, TF, length(grid))
end

export cheb_analysis_matrix
