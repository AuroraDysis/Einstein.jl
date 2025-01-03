@testitem "Continued Fraction" begin
    using GRSuite.QNM, Test

    for type in [Float64, BigFloat]
        max_error = 10 * eps(type)
        tol = eps(type)
        min_iter = 50
        max_iter = 1000
    
        @testset "$type, sqrt(2)" begin
            a(j) = one(type)
            b(j) = j == 0 ? one(type) : type(2)
            f, error, j = contfrac_lentz(type, a, b, tol, min_iter, max_iter)
            @test f ≈ sqrt(type(2))
            @test error <= max_error
        end
    
        @testset "$type, Golden Ratio" begin
            a(j) = one(type)
            b(j) = one(type)
            f, error, j = contfrac_lentz(type, a, b, tol, min_iter, max_iter)
            @test f ≈ (1 + sqrt(type(5))) / 2
            @test error <= max_error
        end
    end
end
