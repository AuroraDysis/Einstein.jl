@testitem "sws_eigvals" begin
    using GRSuite, Test

    # QNM
    a = 0.7
    s = 2
    l = 2
    m = 2
    ω = 0.532600243551018 - 0.08079287315500905im
    eres = -1.096833018126064 + 0.18315825271747369im
    evals = sws_eigvals(Float64, s, a * ω, m, 20)
    eidx = sws_eigvalidx(s, l, m)
    @test evals[eidx] ≈ eres
end

@testitem "SWSFun" begin
    using Integrals, GRSuite, Test

    a = 0.7
    s = 2
    l = 2
    m = 2
    ω = 0.532600243551018 - 0.08079287315500905im

    tol = 10 * eps(Float64)
    f = SWSFun{Float64}(s, a * ω, m, l, 25)
    fi(x, p) = f(x) * conj(f(x)) * sin(x)
    domain = (0, π)
    prob = IntegralProblem(fi, domain)
    sol = solve(prob, QuadGKJL())
    @test isapprox(sol.u, 1, atol=tol)
end
