struct ChebshevInterpolator{TF<:AbstractFloat}
    grid::ChebyshevGrid{TF}
    weights::Vector{TF}  # Barycentric weights
end

function (interpolator::ChebshevInterpolator{TF})(
    values::AbstractVector{TR}, x::TF
) where {TF<:AbstractFloat,TR<:Union{TF,Complex{TF}}}
    return bary(interpolator.weights, interpolator.grid.data, values, x)
end

function cheb_interp(grid::ChebyshevGrid{TF}) where {TF<:AbstractFloat}
    if grid.type == ChebyshevNode.FirstKind
        bary_weights = cheb1_barywts(TF, length(grid))
    elseif grid.type == ChebyshevNode.SecondKind
        bary_weights = cheb2_barywts(TF, length(grid))
    else
        throw(ArgumentError("Invalid Chebyshev node type: $(grid.type)"))
    end

    return ChebshevInterpolator(grid, bary_weights)
end

export ChebshevInterpolator, cheb_interp
