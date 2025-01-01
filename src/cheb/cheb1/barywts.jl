"""
    cheb1_barywts([TR=Float64], n::Integer)

Compute the barycentric weights for Chebyshev points of the first kind.

# Arguments
- `TR`: Type parameter for the weights (e.g., Float64)
- `n`: Number of points

# Returns
- Vector of n barycentric weights

# Mathematical Details
For Chebyshev points of the first kind, the barycentric weights are:
```math
w_j = (-1)^j \\sin\\left(\\frac{(2j+1)\\pi}{2n}\\right), \\quad j = 0,\\ldots,n-1
```

These weights provide optimal stability for barycentric interpolation at
Chebyshev points of the first kind.

See also: [`bary`](@ref), [`cheb1_pts`](@ref)
"""
function cheb1_barywts(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    # Handle corner cases
    if n == 0
        return TR[]
    elseif n == 1
        return [one(TR)]
    end

    half = one(TR) / 2
    pi_over_n = convert(TR, π) / n
    v = zeros(TR, n)
    @inbounds for j in 0:(n - 1)
        θ = (n - j - half) * pi_over_n
        v[j + 1] = sin(θ)
    end

    # The following flipping trick forces symmetry. Also due to the nature of 
    # the sine function, those computed with a big argument are replaced by ones
    # with a small argument, improving the relative accuracy.
    half_n = floor(Int, n / 2)
    # Copy values from end to beginning for symmetry
    @inbounds for i in 1:half_n
        v[i] = v[n - i + 1]
    end

    # Flip signs for odd indices
    @inbounds for i in (n - 1):-2:1
        v[i] = -v[i]
    end

    return v
end

function cheb1_barywts(n::TI) where {TI<:Integer}
    return cheb1_barywts(Float64, n)
end

export cheb1_barywts

@testset "cheb1_barywts" begin
    # Test n=0 case
    @test cheb1_barywts(0) == Float64[]

    # Test n=1 case
    @test cheb1_barywts(1) ≈ [1.0]

    # Test n=5 case
    w5 = cheb1_barywts(5)
    @test w5 ≈ [
        0.309016994374947, -0.809016994374948, 1.0, -0.809016994374948, 0.309016994374947
    ]

    w6 = cheb1_barywts(6)
    @test w6 ≈ [
        -0.258819045102521,
        0.707106781186548,
        -0.965925826289068,
        0.965925826289068,
        -0.707106781186548,
        0.258819045102521,
    ]
end
