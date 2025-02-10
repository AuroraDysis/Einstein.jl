using TestItems

@testitem "cheb1_angles" begin
    using ApproxFun, Einstein.ChebyshevSuite, Test

    for TR in [Float64, BigFloat]
        tol = 10 * eps(TR)
        dom = -one(TR) .. one(TR)
        @testset "$TR" begin
            for n in 0:10
                @test isapprox(
                    cheb1_angles(TR, n), reverse(acos.(points(Chebyshev(dom), n))), atol=tol
                )
            end
        end
    end
end
