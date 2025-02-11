abstract type AbstractChebyshevSynthesisImplementation end

"""
    ChebyshevSynthesis(grid::ChebyshevGrid{TF}) where {TF<:AbstractFloat}

Construct a Chebyshev synthesis operator for the given grid, converting from coefficients to values.
"""
struct ChebyshevSynthesis{
    TF<:AbstractFloat,
    TNode<:AbstractChebyshevNode,
    TSynthesisImpl<:AbstractChebyshevSynthesisImplementation,
}
    grid::ChebyshevGrid{TF,TNode}
    impl::TSynthesisImpl

    function ChebyshevSynthesis(
        grid::ChebyshevGrid{TF,TNode}
    ) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode}
        synthesis = cheb_synthesis(grid.node, TF, length(grid))
        return new{TF,TNode,typeof(synthesis)}(grid, synthesis)
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
