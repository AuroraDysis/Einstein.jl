"""
    ChebyshevInterpolation(grid::ChebyshevGrid{TF}) where {TF<:AbstractFloat}

Construct a barycentric interpolation with precomputed weights for a Chebyshev grid.
"""
struct ChebyshevInterpolation{TF<:AbstractFloat}
    grid::ChebyshevGrid{TF}
    bary_itp::BarycentricInterpolation{TF}

    function ChebyshevInterpolation(grid::ChebyshevGrid{TF}) where {TF<:AbstractFloat}
        if grid.type == ChebyshevNode.FirstKind
            weights = cheb1_barywts(TF, length(grid))
        elseif grid.type == ChebyshevNode.SecondKind
            weights = cheb2_barywts(TF, length(grid))
        else
            throw(ArgumentError("Invalid Chebyshev node type: $(grid.type)"))
        end

        return new{TF}(grid, BarycentricInterpolation(grid.data, weights))
    end
end

function (itp::ChebyshevInterpolation{TF})(
    f::AbstractVector{TR}, x0::TF
) where {TF<:AbstractFloat,TR<:Union{TF,Complex{TF}}}
    return itp.bary_itp(f, x0)
end

export ChebyshevInterpolation
