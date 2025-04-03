@testitem "qnm_pep_companion" begin
    using LinearAlgebra

    # Test 1: Simple quadratic PEP (degree 2) with 2x2 matrices
    A0 = [1.0 2.0; 3.0 4.0]
    A1 = [5.0 6.0; 7.0 8.0]
    A2 = [9.0 10.0; 11.0 12.0]
    pep = [A0, A1, A2]
    
    A, E = qnm_pep_companion(pep)
    
    # Check dimensions
    @test size(A) == (4, 4)  # 2*2 x 2*2 for degree 2
    @test size(E) == (4, 4)
    
    # Check structure of A
    @test A[1:2, 1:2] ≈ -A1
    @test A[1:2, 3:4] ≈ -A0
    @test A[3:4, 1:2] ≈ [1.0 0.0; 0.0 1.0]
    @test A[3:4, 3:4] ≈ [0.0 0.0; 0.0 0.0]
    
    # Check structure of E
    @test E[1:2, 1:2] ≈ A2
    @test E[1:2, 3:4] ≈ [0.0 0.0; 0.0 0.0]
    @test E[3:4, 1:2] ≈ [0.0 0.0; 0.0 0.0]
    @test E[3:4, 3:4] ≈ [1.0 0.0; 0.0 1.0]
    
    # Test 2: Linear PEP (degree 1) with 3x3 matrices
    B0 = rand(3, 3)
    B1 = rand(3, 3)
    pep_linear = [B0, B1]
    
    A_linear, E_linear = qnm_pep_companion(pep_linear)
    
    # Check dimensions
    @test size(A_linear) == (3, 3)  # 3x3 for degree 1
    @test size(E_linear) == (3, 3)
    
    # Check structure of A_linear
    @test A_linear ≈ -B0
    
    # Check structure of E_linear
    @test E_linear ≈ B1
    
    # Test 3: Complex matrices
    C0 = rand(Complex{Float64}, 2, 2)
    C1 = rand(Complex{Float64}, 2, 2)
    C2 = rand(Complex{Float64}, 2, 2)
    pep_complex = [C0, C1, C2]
    
    A_complex, E_complex = qnm_pep_companion(pep_complex)
    
    # Check dimensions
    @test size(A_complex) == (4, 4)
    @test size(E_complex) == (4, 4)
    
    # Check that complex arithmetic worked correctly
    @test A_complex[1:2, 1:2] ≈ -C1
    @test A_complex[1:2, 3:4] ≈ -C0
    @test E_complex[1:2, 1:2] ≈ C2
end

