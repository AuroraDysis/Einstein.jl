abstract type AbstractChebyshevSynthesisImplementation end

struct ChebyshevSynthesis{
    TF<:AbstractFloat,TSynthesisImpl<:AbstractChebyshevSynthesisImplementation
}
    grid::ChebyshevGrid{TF}
    synthesis::TSynthesisImpl

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
