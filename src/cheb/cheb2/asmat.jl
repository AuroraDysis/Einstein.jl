function cheb2_amat(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    A = Array{Float64,2}(undef, n, n)
    op = Cheb2Vals2CoeffsOp(TR, n)
    @inbounds for i in 1:n
        A[:, i] = op(OneElement(one(TR), i, n))
    end
    return A
end

function cheb2_amat(n::TI) where {TI<:Integer}
    return cheb2_amat(Float64, n)
end

function cheb2_smat(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    S = Array{Float64,2}(undef, n, n)
    op = Cheb2Coeffs2ValsOp(TR, n)
    @inbounds for i in 1:n
        S[:, i] = op(OneElement(one(TR), i, n))
    end
    return S
end

function cheb2_smat(n::TI) where {TI<:Integer}
    return cheb2_smat(Float64, n)
end

export cheb2_amat, cheb2_smat

@testset "cheb2_amat, cheb2_smat" begin
    n = 32  # Enough points for good accuracy
    x = cheb2_pts(Float64, n)
    A = cheb2_amat(Float64, n)
    S = cheb2_smat(Float64, n)

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
