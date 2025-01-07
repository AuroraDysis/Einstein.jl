@testitem "qnm_kerr_radial_di" begin
    using StaticArrays, GRSuite, Test

    a = 0.7
    s = 2
    l = 2
    m = 2
    ω = 0.532600243551018 - 0.08079287315500905im
    A = -1.096833018126064 + 0.18315825271747035im
    Λ = -A
    tol = 1e-12

    params = QNMKerrRadialDIParams{Float64}(; a=a, s=s, m=m, ω_guess=ω)
    cache = QNMKerrRadialDICache{Float64}(params, ω, Λ)
    δ = qnm_kerr_radial_di_δ(SA[real(ω), imag(ω)], cache)

    @test abs(δ[1]) < tol
    @test abs(δ[2]) < tol
end

@testitem "qnm_kerr_cf" begin
    a = 0.7
    s = 2
    l = 2
    m = 2
    n = 0
    ω = 0.532600243551018 - 0.08079287315500766im
    A = -1.096833018126064 + 0.18315825271747035im
    l_max = 20

    ω_pert = ω + rand(Complex{Float64}) / 1000

    params = QNMKerrCFParams{Float64}(; a=a, s=s, l=l, m=m, n=n, ω_guess=ω_pert, l_max=l_max)

    ωsol = qnm_kerr_cf(params)

    tol = typetol(Float64)
    @test abs(ωsol - ω) < tol
end

@testitem "qnm_kerr" begin
    a = 0.7
    s = 2
    l = 2
    m = 2
    n = 0
    ω = 0.532600243551018 - 0.08079287315500766im
    A = -1.096833018126064 + 0.18315825271747035im
    l_max = 20

    ω_pert = ω + rand(Complex{Float64}) / 1000

    params = QNMKerrChebParams{Float64}(;
        a=a, s=s, l=l, m=m, n=n, ω_guess=ω_pert, l_max=l_max, cheb_n=80
    )

    ωsol = qnm_kerr_cheb(params)

    tol = typetol(Float64)
    @test abs(ωsol - ω) < tol
end
