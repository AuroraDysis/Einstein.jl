abstract type AbstractChebyshevSynthesisImplementation end

"""
    ChebyshevSynthesis(grid::ChebyshevGrid)

Construct a Chebyshev synthesis operator for the given grid, converting from coefficients to values.
"""
struct ChebyshevSynthesis{
    TF<:AbstractFloat,
    Basis<:AbstractBasisFunction,
    TSynthesisImpl<:AbstractChebyshevSynthesisImplementation,
}
    grid::ChebyshevGrid{TF,Basis}
    impl::TSynthesisImpl

    function ChebyshevSynthesis(
        grid::ChebyshevGrid{TF,Basis}
    ) where {TF<:AbstractFloat,Basis<:AbstractBasisFunction}
        impl = _cheb_synthesis(grid.basis, TF, length(grid))
        return new{TF,Basis,typeof(impl)}(grid, impl)
    end
end

function (syn::ChebyshevSynthesis{TF,Basis,TSynthesisImpl})(
    coeffs::AbstractVector{TFC}
) where {
    TF<:AbstractFloat,
    Basis<:AbstractBasisFunction,
    TSynthesisImpl<:AbstractChebyshevSynthesisImplementation,
    TFC<:Union{TF,Complex{TF}},
}
    return syn.impl(coeffs)
end

export ChebyshevSynthesis
