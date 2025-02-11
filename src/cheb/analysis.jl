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
        analysis = cheb_analysis(grid.node, TF, length(grid))
        return new{TF,TNode,typeof(analysis)}(grid, analysis)
    end
end

function (syn::ChebyshevAnalysis{TF,TAnalysisImpl})(
    coeffs::AbstractVector{TRC}
) where {
    TF<:AbstractFloat,
    TRC<:Union{TF,Complex{TF}},
    TAnalysisImpl<:AbstractChebyshevAnalysisImplementation,
}
    return impl.analysis(coeffs)
end

export ChebyshevAnalysis
