using TestItems

@testitem "chebtech1_angles" begin
    using ApproxFun, Einstein.ChebyshevSuite, Test

    for TF in [Float64, BigFloat]
        tol = 10 * eps(TF)
        dom = -one(TF) .. one(TF)
        @testset "$TF" begin
            for n in 0:10
                @test isapprox(
                    chebtech1_angles(TF, n), reverse(acos.(points(Chebyshev(dom), n))), atol=tol
                )
            end
        end
    end
end
