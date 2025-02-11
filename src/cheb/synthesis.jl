abstract type AbstractChebyshevSynthesisImplementation end

"""
    ChebyshevSynthesis(grid::ChebyshevGrid)

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
        impl = _cheb_synthesis(grid.node, TF, length(grid))
        return new{TF,TNode,typeof(impl)}(grid, impl)
    end
end

function (syn::ChebyshevSynthesis{TF,TNode,TSynthesisImpl})(
    coeffs::AbstractVector{TFC}
) where {
    TF<:AbstractFloat,
    TNode<:AbstractChebyshevNode,
    TSynthesisImpl<:AbstractChebyshevSynthesisImplementation,
    TFC<:Union{TF,Complex{TF}},
}
    return syn.impl(coeffs)
end

export ChebyshevSynthesis
