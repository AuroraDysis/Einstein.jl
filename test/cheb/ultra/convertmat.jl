using TestItems

@testitem "ultra_convertmat" begin
    using ApproxFun, PDESuite.ChebSuite, Test

    for type in [Float64, BigFloat]
        tol = 10 * eps(type)
        for n in 2:10
            for K1 in 1:3
                for K2 in K1:3
                    dom = -one(type) .. one(type)
                    @test isapprox(
                        ultra_convertmat(type, n, K1, K2),
                        Conversion(Ultraspherical(K1, dom), Ultraspherical(K2, dom))[
                            1:n, 1:n
                        ],
                        atol=tol,
                    )
                end
            end
        end
    end
end
