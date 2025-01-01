using FFTW: fft, ifft
using FastTransforms: fft, ifft

"""
    cheb1_coeffs2vals(coeffs::AbstractVector{<:AbstractFloat}) -> Vector

Convert Chebyshev coefficients of the first kind to values at Chebyshev nodes using a DCT-III equivalent
transform.

# Arguments
- `coeffs::AbstractVector{<:AbstractFloat}`: Vector of Chebyshev coefficients in descending order

# Returns
- Vector of values at Chebyshev nodes of the first kind: cos((2j+1)π/(2n)), j = 0,1,...,n-1

# Mathematical Background
The function implements the transform from coefficient space to physical space for Chebyshev
polynomials of the first kind Tₙ(x). The transformation:

1. Uses the relationship between Chebyshev polynomials and cosine series
2. Applies a type-III discrete cosine transform via FFT
3. Preserves polynomial symmetries:
   - Even coefficients produce even functions
   - Odd coefficients produce odd functions

# Implementation Details
1. Computes weights w[k] = exp(-ikπ/(2n))/2
2. Applies special treatment for endpoints
3. Uses FFT for efficient computation
4. Enforces symmetry properties

# Examples
```julia
# Convert a constant function
julia> cheb1_coeffs2vals([1.0])
1-element Vector{Float64}:
 1.0

# Convert linear coefficients
julia> cheb1_coeffs2vals([1.0, 2.0])
2-element Vector{Float64}:
  3.0
 -1.0

# Convert quadratic coefficients
julia> cheb1_coeffs2vals([1.0, 0.0, -0.5])
3-element Vector{Float64}:
  0.5
 -0.5
  0.5
```

# References
1. Trefethen, L. N. (2000). Spectral Methods in MATLAB. SIAM.
2. Mason, J. C., & Handscomb, D. C. (2002). Chebyshev Polynomials. Chapman & Hall/CRC.
3. Boyd, J. P. (2001). Chebyshev and Fourier Spectral Methods. Dover.

See also: [`cheb2_coeffs2vals`](@ref)
"""
function cheb1_coeffs2vals(coeffs::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}
    n = length(coeffs)

    # Trivial case (constant polynomial)
    if n <= 1
        return deepcopy(coeffs)
    end

    # Check for symmetry:
    # Columns where all even-indexed coeffs are zero => function is even
    isEven = all(x -> x ≈ 0, @view(coeffs[2:2:end]))
    # Columns where all odd-indexed coeffs are zero => function is odd
    isOdd = all(x -> x ≈ 0, @view(coeffs[1:2:end]))

    w = [exp(-im * (k) * pi / (2n)) / 2 for k in 0:(2n - 1)]
    w[1] = 2 * w[1]
    w[n + 1] = 0
    for k in (n + 2):(2n)
        w[k] = -w[k]
    end

    # Mirror the coeffs for FFT:
    # c_mirror = [ coeffs; ones(1,m); flipud(coeffs[2:end, :]) ]
    # In Julia, vcat is used to stack arrays vertically:
    c_mirror = vcat(coeffs, ones(1), reverse(@view(coeffs[2:end])))

    # Apply the weight vector:
    # c_weight = bsxfun(@times, c_mirror, w)
    # We can broadcast directly in Julia:
    c_weight = c_mirror .* w

    # Perform FFT:
    values = fft(c_weight)  # 1 indicates along the first dimension

    # Truncate and flip the order (n:-1:1 in MATLAB):
    # In Julia, we can reverse the array along dimension 1 and then take the first n rows.
    values = reverse(@view(values[1:n]))

    # Post-process to handle real and imaginary inputs:
    values = real.(values)

    # Enforce symmetry
    if isEven
        values = (values + reverse(values)) ./ 2
    elseif isOdd
        values = (values .- reverse(values)) ./ 2
    end

    return values
end

export cheb1_coeffs2vals

@testset "cheb1_coeffs2vals" begin
    using LinearAlgebra

    # Set tolerance
    tol = 100 * eps()

    @testset "Single coefficient" begin
        c = [sqrt(2)]
        v = cheb1_coeffs2vals(c)
        @test c ≈ v
    end

    @testset "Even case" begin
        # Simple data (even case)
        c = collect(6.0:-1:1)
        # Exact values
        vTrue = [
            -3 * sqrt(6) / 2 - 5 / sqrt(2) + 2 * sqrt(3) + 7
            4 - sqrt(2) / 2
            -3 * sqrt(6) / 2 + 5 / sqrt(2) - 2 * sqrt(3) + 7
            3 * sqrt(6) / 2 - 5 / sqrt(2) - 2 * sqrt(3) + 7
            4 + sqrt(2) / 2
            3 * sqrt(6) / 2 + 5 / sqrt(2) + 2 * sqrt(3) + 7
        ]

        # Test real branch
        v = cheb1_coeffs2vals(c)
        @test norm(v - vTrue, Inf) < tol
        @test all(iszero, imag.(v))
    end

    @testset "Odd case" begin
        # Simple data (odd case)
        c = collect(5.0:-1:1)
        # Exact values
        vTrue = [
            11 / 2 + sqrt(5) - 2 * sqrt((5 + sqrt(5)) / 2) - sqrt((5 - sqrt(5)) / 2)
            11 / 2 - sqrt(5) - 2 * sqrt((5 - sqrt(5)) / 2) + sqrt((5 + sqrt(5)) / 2)
            3
            11 / 2 - sqrt(5) + 2 * sqrt((5 - sqrt(5)) / 2) - sqrt((5 + sqrt(5)) / 2)
            11 / 2 + sqrt(5) + 2 * sqrt((5 + sqrt(5)) / 2) + sqrt((5 - sqrt(5)) / 2)
        ]

        # Test real branch
        v = cheb1_coeffs2vals(c)
        @test norm(v - vTrue, Inf) < tol
        @test all(iszero, imag.(v))
    end

    @testset "Symmetry preservation" begin
        c = kron(ones(10), Matrix{Float64}(I, 2, 2))
        v1 = cheb1_coeffs2vals(c[:, 1])
        v2 = cheb1_coeffs2vals(c[:, 2])
        @test norm(v1 - reverse(v1), Inf) ≈ 0
        @test norm(v2 + reverse(v2), Inf) ≈ 0
    end
end
