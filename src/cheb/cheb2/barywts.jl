"""
    cheb2_barywts([TR=Float64], n::Integer)

Compute the barycentric weights for Chebyshev points of the second kind.

# Arguments
- `TR`: Type parameter for the weights (e.g., Float64)
- `n`: Number of points

# Returns
- Vector of n barycentric weights

# Mathematical Details
For Chebyshev points of the second kind, the barycentric weights are:
```math
w_j = (-1)^j \\delta_j, \\quad j = 0,\\ldots,n-1
```
where ``\\delta_j`` is defined as:
```math
\\delta_j = \\begin{cases}
1/2 & j = 0 \\text{ or } j = n-1 \\\\
1 & \\text{otherwise}
\\end{cases}
```

These weights are optimized for numerical stability and efficiency in the barycentric
interpolation formula.

See also: [`bary`](@ref), [`cheb2_pts`](@ref)
"""
function cheb2_barywts(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    if n == 0
        return TR[]
    elseif n == 1
        return [one(TR)]
    end

    w = ones(TR, n)

    @inbounds begin
        w[(end - 1):-2:1] .= -1
        w[1] /= 2
        w[end] = one(TR) / 2
    end

    return w
end

function cheb2_barywts(n::TI) where {TI<:Integer}
    return cheb2_barywts(Float64, n)
end

export cheb2_barywts

@testset "cheb2_barywts" begin
    # Test n=1 case
    w1 = cheb2_barywts(Float64, 1)
    @test w1 ≈ [1.0]

    # Test n=2 case
    w2 = cheb2_barywts(Float64, 2)
    @test w2 ≈ [-0.500000000000000, 0.500000000000000]

    # Test n=5 case
    w5 = cheb2_barywts(Float64, 5)
    @test w5 ≈ [0.500000000000000, -1.0, 1.0, -1.0, 0.500000000000000]

    w6 = cheb2_barywts(Float64, 6)
    @test w6 ≈ [-0.500000000000000, 1.0, -1.0, 1.0, -1.0, 0.500000000000000]
end
