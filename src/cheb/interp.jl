struct ChebshevInterpolator{TF<:AbstractFloat,TNode<:AbstractChebyshevNode}
    grid::ChebyshevGrid{TF,TNode}
    weights::Vector{TF}  # Barycentric weights
end

function (interpolator::ChebshevInterpolator{TF,TNode})(
    values::AbstractVector{TR}, x::TF
) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode,TR<:Union{TF,Complex{TF}}}
    return bary(interpolator.weights, interpolator.grid.data, values, x)
end

function cheb_interp(
    grid::ChebyshevGrid{TF,TNode}
) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode}
    if TNode == ChebyshevNode.FirstKind
        bary_weights = cheb1_barywts(TF, length(grid))
    elseif TNode == ChebyshevNode.SecondKind
        bary_weights = cheb2_barywts(TF, length(grid))
    else
        throw(ArgumentError("Invalid Chebyshev node type: $TNode"))
    end

    return ChebshevInterpolator(grid, bary_weights)
end

export ChebshevInterpolator, cheb_interp
