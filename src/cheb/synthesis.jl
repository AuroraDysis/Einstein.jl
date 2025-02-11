abstract type AbstractChebyshevSynthesisImplementation end

"""
    ChebyshevSynthesis(grid::ChebyshevGrid{TF}) where {TF<:AbstractFloat}

Construct a Chebyshev synthesis operator for the given grid, converting from coefficients to values.
"""
struct ChebyshevSynthesis{
    TF<:AbstractFloat,TSynthesisImpl<:AbstractChebyshevSynthesisImplementation
}
    grid::ChebyshevGrid{TF}
    impl::TSynthesisImpl

    function ChebyshevSynthesis(grid::ChebyshevGrid{TF}) where {TF<:AbstractFloat}
        if grid.type == ChebyshevNode.FirstKind
            synthesis = ChebyshevFirstKindSynthesis{TF}(length(grid))
        elseif grid.type == ChebyshevNode.SecondKind
            synthesis = ChebyshevSecondKindSynthesis{TF}(length(grid))
        else
            throw(
                ArgumentError(
                    "kind must be either ChebyshevNode.FirstKind or ChebyshevNode.SecondKind",
                ),
            )
        end

        return new{TF,typeof(synthesis)}(grid, synthesis)
    end
end

function (syn::ChebyshevSynthesis{TF,TSynthesisImpl})(
    coeffs::AbstractVector{TRC}
) where {
    TF<:AbstractFloat,
    TRC<:Union{TF,Complex{TF}},
    TSynthesisImpl<:AbstractChebyshevSynthesisImplementation,
}
    return impl.synthesis(coeffs)
end

export ChebyshevSynthesis
