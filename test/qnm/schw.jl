@testitem "qnm_schw_chebnep" begin
    using LinearAlgebra, GRSuite, Test

    s = ℓ = 2
    n = 80
    ω = 0.3736716844180419 - 0.0889623156889357im
    tol = 1e-14

    @testset "ReggeWheeler" begin
        pep = qnm_schw_chebnep(Float64, s, ℓ, n, SchwPType.ReggeWheeler)
        A, E = qnm_polyeig(ComplexF64, pep)
        # Ax = λEx
        λ = eigvals!(A, E; sortby=abs)
        λclosest = argmin(x -> abs(x - ω), λ)
        @test isapprox(λclosest, ω, atol=tol)
    end

    @testset "Zerilli" begin
        pep = qnm_schw_chebnep(Float64, s, ℓ, n, SchwPType.Zerilli)
        A, E = qnm_polyeig(ComplexF64, pep)
        λ = eigvals!(A, E; sortby=abs)
        λclosest = argmin(x -> abs(x - ω), λ)
        @test isapprox(λclosest, ω, atol=tol)
    end
end
