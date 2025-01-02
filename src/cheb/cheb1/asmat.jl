"""
    cheb1_amat([TR=Float64], n::Integer) -> Matrix{TR}

Construct the analysis matrix A that transforms function values at Chebyshev points of the first kind to Chebyshev coefficients.

# Arguments
- `TR`: Element type (defaults to Float64)
- `n`: Number of points/coefficients
"""
function cheb1_amat(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    A = Array{Float64,2}(undef, n, n)
    op = Cheb1Vals2CoeffsOp(TR, n)
    @inbounds for i in 1:n
        A[:, i] = op(OneElement(one(TR), i, n))
    end
    return A
end

"""
    cheb1_smat([TR=Float64], n::Integer) -> Matrix{TR}

Construct the synthesis matrix S that transforms Chebyshev coefficients to function values at Chebyshev points of the first kind.

# Arguments
- `TR`: Element type (defaults to Float64)
- `n`: Number of points/coefficients
"""
function cheb1_smat(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    S = Array{Float64,2}(undef, n, n)
    op = Cheb1Coeffs2ValsOp(TR, n)
    @inbounds for i in 1:n
        S[:, i] = op(OneElement(one(TR), i, n))
    end
    return S
end

function cheb1_amat(n::TI) where {TI<:Integer}
    return cheb1_amat(Float64, n)
end

function cheb1_smat(n::TI) where {TI<:Integer}
    return cheb1_smat(Float64, n)
end

export cheb1_amat, cheb1_smat

@testset "cheb1_amat, cheb1_smat" begin
    n = 32  # Enough points for good accuracy
    x = cheb1_pts(Float64, n)
    A = cheb1_amat(Float64, n)
    S = cheb1_smat(Float64, n)

    @testset "Transform and recover" begin
        # Test with polynomial that should be exactly represented
        f = @. 3x^2 + 2x - 1
        coeffs = A * f
        f_recovered = S * coeffs
        @test f_recovered ≈ f rtol = 1e-12

        # Test with trigonometric function
        f = @. sin(π * x) * cos(2π * x)
        coeffs = A * f
        f_recovered = S * coeffs
        @test f_recovered ≈ f rtol = 1e-12
    end
end
