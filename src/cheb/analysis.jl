abstract type AbstractChebyshevAnalysisImplementation end

"""
    ChebyshevAnalysis(grid::ChebyshevGrid)

Construct a Chebyshev analysis operator for the given grid, converting from coefficients to values.
"""
struct ChebyshevAnalysis{
    TF<:AbstractFloat,
    TNode<:AbstractChebyshevNode,
    TAnalysisImpl<:AbstractChebyshevAnalysisImplementation,
}
    grid::ChebyshevGrid{TF,TNode}
    impl::TAnalysisImpl

    function ChebyshevAnalysis(
        grid::ChebyshevGrid{TF,TNode}
    ) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode}
        impl = cheb_analysis(grid.node, TF, length(grid))
        return new{TF,TNode,typeof(impl)}(grid, impl)
    end
end

function (ana::ChebyshevAnalysis{TF,TNode,TAnalysisImpl})(
    values::AbstractVector{TFC}
) where {
    TF<:AbstractFloat,
    TNode<:AbstractChebyshevNode,
    TAnalysisImpl<:AbstractChebyshevAnalysisImplementation,
    TFC<:Union{TF,Complex{TF}},
}
    return ana.impl(values)
end

export ChebyshevAnalysis
