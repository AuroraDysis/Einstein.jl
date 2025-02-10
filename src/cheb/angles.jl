@doc raw"""
    cheb_angles(grid::ChebyshevGrid{TF}) where {TF<:AbstractFloat}

Compute angles for Chebyshev grid of given type.
"""
function cheb_angles(grid::ChebyshevGrid{TF}) where {TF<:AbstractFloat}
    if grid.type == ChebyshevNode.FirstKind
        return cheb1_angles(TF, length(grid))
    elseif grid.type == ChebyshevNode.SecondKind
        return cheb2_angles(TF, length(grid))
    else
        throw(ArgumentError("Invalid Chebyshev node type: $(grid.type)"))
    end
end

export cheb_angles
