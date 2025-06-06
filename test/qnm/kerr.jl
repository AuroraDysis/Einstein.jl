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

    ωsol = qnm_kerr(a, s, l, m, ω_pert; n=n, l_max=l_max)

    tol = typetol(Float64)
    @test abs(ωsol - ω) < tol
end
