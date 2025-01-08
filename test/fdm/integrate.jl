@testitem "fdm_integrate_simpson" begin
    using GRSuite.FDMSuite, Test

    dx = 0.01
    x = fdm_grid(0.0, 1.0, dx)
    fx = @. sin(2π * x)
    integral = fdm_integrate_simpson(fx, dx)
    @test integral ≈ 0.0 atol = 1e-14
end

@testitem "fdm_integrate_trapezoidal" begin
    using GRSuite.FDMSuite, Test

    dx = 0.01
    x = fdm_grid(0.0, 1.0, dx)
    fx = @. sin(2π * x)
    integral = fdm_integrate_trapezoidal(fx, dx)
    @test integral ≈ 0.0 atol = 1e-14
end
