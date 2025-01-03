@testset "ultra_multmat" begin
    for TR in [Float64, BigFloat]
        tol = 1000 * eps(TR)
        for n in 5:10
            dom = -one(TR) .. one(TR)
            coeffs = collect(TR, 1:n)
            f = Fun(Chebyshev(dom), coeffs)
            for λ in 0:3
                t1 = ultra_multmat(coeffs, λ)
                t2 = if λ == 0
                    Matrix(@view(Multiplication(f, Chebyshev(dom))[1:n, 1:n]))
                else
                    Matrix(@view(Multiplication(f, Ultraspherical(λ, dom))[1:n, 1:n]))
                end
                @test isapprox(t1, t2, atol=tol)
            end
        end
    end
end
