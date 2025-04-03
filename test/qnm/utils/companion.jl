@testitem "qnm_pep_companion" begin
    using LinearAlgebra, SparseArrays

    @testset "Simple quadratic PEP (degree 2) with 2x2 matrices" begin
        A0 = [1 2; 3 4]
        A1 = [5 6; 7 8]
        A2 = [9 10; 11 12]
        pep = [A0, A1, A2]

        A_exact = [-5 -6 -1 -2; -7 -8 -3 -4; 1 0 0 0; 0 1 0 0]
        E_exact = [9 10 0 0; 11 12 0 0; 0 0 1 0; 0 0 0 1]

        A, E = qnm_pep_companion(pep)

        @test size(A) == (4, 4)
        @test size(E) == (4, 4)

        @test isapprox(A, A_exact)
        @test isapprox(E, E_exact)
    end

    @testset "Linear PEP (degree 1) with 3x3 matrices" begin
        B0 = rand(3, 3)
        B1 = rand(3, 3)
        pep_linear = [B0, B1]

        A_linear, E_linear = qnm_pep_companion(pep_linear)

        @test size(A_linear) == (3, 3)  # 3x3 for degree 1
        @test size(E_linear) == (3, 3)

        @test A_linear ≈ -B0

        @test E_linear ≈ B1
    end

    @testset "Complex matrices" begin
        A0 = [1+1im 2+2im; 3+3im 4+4im]
        A1 = [5+5im 6+6im; 7+7im 8+8im]
        A2 = [9+9im 10+10im; 11+11im 12+12im]
        pep_complex = [A0, A1, A2]

        A_complex_exact = [
            -5-5im -6-6im -1-1im -2-2im
            -7-7im -8-8im -3-3im -4-4im
            1+0im 0+0im 0+0im 0+0im
            0+0im 1+0im 0+0im 0+0im
        ]
        E_complex_exact = [
            9+9im 10+10im 0+0im 0+0im
            11+11im 12+12im 0+0im 0+0im
            0+0im 0+0im 1+0im 0+0im
            0+0im 0+0im 0+0im 1+0im
        ]

        A_complex, E_complex = qnm_pep_companion(pep_complex)

        @test size(A_complex) == (4, 4)
        @test size(E_complex) == (4, 4)

        @test isapprox(A_complex, A_complex_exact)
        @test isapprox(E_complex, E_complex_exact)
    end

    @testset "Sparse matrices" begin
        S0 = sprand(ComplexF64, 10, 10, 0.1)
        S1 = sprand(ComplexF64, 10, 10, 0.1)
        S2 = sprand(ComplexF64, 10, 10, 0.1)
        pep_sparse = [S0, S1, S2]

        A_sparse, E_sparse = qnm_pep_companion(pep_sparse)

        @test issparse(A_sparse)
        @test issparse(E_sparse)
    end
end
