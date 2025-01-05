@testitem "qnm_kerrnep" begin
    using GRSuite, Test

    a = 0.7
    s = 2
    l = 2
    m = 2
    n = 80
    ω = 0.532600243551018 - 0.08079287315500905im
    tol = 1e-12

    cache = qnm_kerrnep_cache(Float64, a, s, m, n)
    δ = qnm_kerrnep_step!(cache, ω, l)

    @test abs(δ) < tol
end
