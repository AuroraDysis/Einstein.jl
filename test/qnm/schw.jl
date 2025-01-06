@testitem "qnm_schw_chebnep" begin
    using GRSuite, Test, NonlinearEigenproblems

    s = ℓ = 2
    n = 80
    ω = 0.3736716844180419 - 0.0889623156889357im
    tol = 1e-14

    @testset "ReggeWheeler" begin
        nep = qnm_schw_chebnep(Float64, s, ℓ, n, SchwPType.ReggeWheeler)
        λ, v = polyeig(nep)
        @test isapprox(λ[1], ω, atol=tol)
    end

    @testset "Zerilli" begin
        nep = qnm_schw_chebnep(Float64, s, ℓ, n, SchwPType.Zerilli)
        λ, v = polyeig(nep)
        @test isapprox(λ[1], ω, atol=tol)
    end
end
