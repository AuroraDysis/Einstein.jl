@doc raw"""
    cheb_angles(grid::ChebyshevGrid{TF,Basis}) where {TF<:AbstractFloat,Basis<:AbstractBasisFunction}

Compute angles for Chebyshev grid of given type.
"""
function cheb_angles(
    grid::ChebyshevGrid{TF,Basis}
) where {TF<:AbstractFloat,Basis<:AbstractBasisFunction}
    return _cheb_angles(grid.basis, TF, length(grid))
end

export cheb_angles
