@testitem "cheb_series_chop" begin
    using LinearAlgebra

    # Basic test cases based on standardChop examples
    
    # Case 1: Simple decreasing coefficients
    coeffs = 10.0 .^ (-(1:50))
    @test cheb_series_chop(coeffs) == 18
    
    # Case 2: Small random perturbation (1e-16)
    random = cos.((1:50) .^ 2)
    perturbed_coeffs = coeffs .+ 1e-16 .* random
    @test cheb_series_chop(perturbed_coeffs) == 15
    
    # Case 3: Medium random perturbation (1e-13)
    perturbed_coeffs = coeffs .+ 1e-13 .* random
    @test cheb_series_chop(perturbed_coeffs) == 13
    
    # Case 4: Large random perturbation (1e-10)
    perturbed_coeffs = coeffs .+ 1e-10 .* random
    @test cheb_series_chop(perturbed_coeffs) == 50
    
    # Case 5: Custom tolerance
    @test cheb_series_chop(coeffs .+ 1e-10 .* random, 1e-10) == 10
    
    # Additional test cases
    
    # Case 6: All zeros should return 1
    @test cheb_series_chop(zeros(20)) == 1
    
    # Case 7: Small vector (< 17 elements) triggering the tol-based implementation
    small_coeffs = 10.0 .^ (-(1:10))
    cutoff = cheb_series_chop(small_coeffs)
    @test 1 ≤ cutoff ≤ 10
    
    # Case 8: Test with complex coefficients
    complex_coeffs = (1.0 .+ 1.0im) .* coeffs
    @test cheb_series_chop(complex_coeffs) == 18
    
    # Case 9: Test with pre-allocated envelope
    envelope = similar(coeffs)
    @test cheb_series_chop(coeffs; envelope=envelope) == 18
end
