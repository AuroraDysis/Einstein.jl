using TestItems

@testitem "ultra_diffmat" begin
    using ApproxFun, Einstein.ChebSuite, Test

    for n in 2:10
        for m in 1:3
            @test isapprox(ultra_diffmat(m, n), Derivative(Chebyshev(), m)[1:n, 1:n])
        end
    end
end
