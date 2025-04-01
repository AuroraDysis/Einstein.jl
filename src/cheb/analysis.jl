abstract type AbstractChebyshevAnalysisImplementation end

"""
    ChebyshevAnalysis(grid::ChebyshevGrid)

Construct a Chebyshev analysis operator for the given grid, converting from coefficients to values.
"""
struct ChebyshevAnalysis{
    TF<:AbstractFloat,
    Basis<:AbstractBasisFunction,
    TAnalysisImpl<:AbstractChebyshevAnalysisImplementation,
}
    grid::ChebyshevGrid{TF,Basis}
    impl::TAnalysisImpl

    function ChebyshevAnalysis(
        grid::ChebyshevGrid{TF,Basis}
    ) where {TF<:AbstractFloat,Basis<:AbstractBasisFunction}
        impl = _cheb_analysis(grid.basis, TF, length(grid))
        return new{TF,Basis,typeof(impl)}(grid, impl)
    end
end

function (ana::ChebyshevAnalysis{TF,Basis,TAnalysisImpl})(
    values::AbstractVector{TFC}
) where {
    TF<:AbstractFloat,
    Basis<:AbstractBasisFunction,
    TAnalysisImpl<:AbstractChebyshevAnalysisImplementation,
    TFC<:Union{TF,Complex{TF}},
}
    return ana.impl(values)
end

export ChebyshevAnalysis
