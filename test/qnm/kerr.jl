@testitem "qnm_kerr_cf" begin
    using GRSuite, Test

    a = 0.7
    s = 2
    l = 2
    m = 2
    n = 0
    ω = 0.532600243551018 - 0.08079287315500766im
    A = -1.096833018126064 + 0.18315825271747035im
    l_max = 20

    ω_pert = ω + rand(Complex{Float64}) / 1000

    params = QNMKerrCFParams{Float64}(;
        a=a, s=s, l=l, m=m, n=n, ω_guess=ω_pert, l_max=l_max
    )

    ωsol = qnm_kerr_cf(params)

    tol = typetol(Float64)
    @test abs(ωsol - ω) < tol
end

@testitem "qnm_kerr_cheb" begin
    using GRSuite, Test

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

@testitem "qnm_kerr_di" begin
    using GRSuite, Test

    a = 0.7
    s = 2
    l = 2
    m = 2
    n = 0
    ω = 0.532600243551018 - 0.08079287315500766im
    A = -1.096833018126064 + 0.18315825271747035im
    l_max = 20

    ω_pert = ω + rand(Complex{Float64}) / 1000

    params = QNMKerrDIParams{Float64}(;
        a=a, s=s, l=l, m=m, n=n, ω_guess=ω_pert, l_max=l_max
    )

    ωsol = qnm_kerr_di(params)

    tol = typetol(Float64)
    @test abs(ωsol - ω) < tol
end
