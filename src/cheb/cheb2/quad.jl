"""
    cheb2_quadwts([TR=Float64], n::TI) where {TR<:AbstractFloat,TI<:Integer}

Compute quadrature weights for Chebyshev points of the 2nd kind.

# Arguments
- `TR`: Type parameter for the weights (e.g., Float64)
- `n`: Number of points

# References
- [chebfun/@chebtech2/quadwts.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/quadwts.m)
"""
function cheb2_quadwts(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    if n == 0
        return TR[]
    elseif n == 1
        return TR[2]
    end

    nm1 = n - 1

    # Fill exact integrals of T_k (even k)
    c = Array{Complex{TR}}(undef, nm1)
    @inbounds begin
        c[1] = TR(2)  # k = 0 case
        for k in 2:2:nm1
            c[k ÷ 2 + 1] = 2 / (one(TR) - k^2)
        end

        # Mirror for DCT via FFT (in-place)
        half_idx = floor(Int, n / 2)
        for i in 2:half_idx
            c[n - i + 1] = c[i]
        end

        # Compute weights via inverse FFT (in-place)
        ifft!(c)
    end

    # Adjust boundary weights (in-place)
    w = Vector{TR}(undef, n)
    @inbounds begin
        w[1] = real(c[1]) / 2
        for i in 2:nm1
            w[i] = real(c[i])
        end
        w[n] = real(c[1]) / 2
    end

    return w
end

function cheb2_quadwts(n::TI) where {TI<:Integer}
    return cheb2_quadwts(Float64, n)
end

export cheb2_quadwts

@testset "cheb2_quadwts" begin
    # Test n=0 case
    @test cheb2_quadwts(0) == Float64[]

    # Test n=1 case
    @test cheb2_quadwts(1) ≈ [2.0]

    # Test n=5 case
    w5 = cheb2_quadwts(5)
    @test w5 ≈ [
        0.0666666666666667,
        0.533333333333333,
        0.800000000000000,
        0.533333333333333,
        0.0666666666666667,
    ]

    w6 = cheb2_quadwts(6)
    @test w6 ≈ [
        0.0400000000000000,
        0.360743041200011,
        0.599256958799989,
        0.599256958799989,
        0.360743041200011,
        0.0400000000000000,
    ]
end
