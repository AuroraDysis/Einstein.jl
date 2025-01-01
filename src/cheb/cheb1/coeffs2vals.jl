"""
    cheb1_coeffs2vals(coeffs::AbstractVector{<:AbstractFloat}) -> Vector
    cheb1_coeffs2vals!(coeffs, cache::Cheb1Coeffs2ValsCache) -> Vector

Convert Chebyshev coefficients of the first kind to values at Chebyshev nodes.

# Performance Guide
For best performance, especially in loops or repeated calls:
1. Create a cache: `cache = Cheb1Coeffs2ValsCache{Float64}(n)`
2. Use the in-place version: `cheb1_coeffs2vals!(coeffs, cache)`

# Arguments
- `coeffs::AbstractVector{<:AbstractFloat}`: Chebyshev coefficients in descending order
- `cache::Cheb1Coeffs2ValsCache`: Pre-allocated workspace (for in-place version)

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
    if n <= 1
        return deepcopy(coeffs)
    end

    cache = Cheb1Coeffs2ValsCache{TR}(n)
    cheb1_coeffs2vals!(coeffs, cache)

    return cache.vals
end

struct Cheb1Coeffs2ValsCache{TR}
    w::Vector{Complex{TR}}    # Weight vector
    tmp::Vector{Complex{TR}}  # Temporary storage
    vals::Vector{TR}         # Result storage

    function Cheb1Coeffs2ValsCache{TR}(n::Integer) where {TR<:AbstractFloat}
        # Precompute weights
        w = Vector{Complex{TR}}(undef, 2n)
        @inbounds begin
            m_im_pi_over_2n = -im * convert(TR, π) / (2n)
            for k in 0:(n - 1)
                w[k + 1] = exp(k * m_im_pi_over_2n) / 2
            end
            w[1] *= 2
            w[n + 1] = 0
            for k in (n + 1):(2n - 1)
                w[k + 1] = -exp(k * m_im_pi_over_2n) / 2
            end
        end
        return new(w, Array{Complex{TR}}(undef, 2n), Vector{TR}(undef, n))
    end
end

function cheb1_coeffs2vals!(
    coeffs::VT, cache::Cheb1Coeffs2ValsCache{TR}
) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}
    n = length(coeffs)
    if n <= 1
        cache.vals .= coeffs
        return cache.vals
    end

    w = cache.w
    tmp = cache.tmp
    vals = cache.vals

    # Check for symmetry
    isEven = all(x -> x ≈ 0, @view(coeffs[2:2:end]))
    isOdd = all(x -> x ≈ 0, @view(coeffs[1:2:end]))

    # Copy coefficients and mirror
    @inbounds begin
        # First half: original coefficients
        for i in 1:n
            tmp[i] = coeffs[i]
        end
        # Second half: mirrored coefficients
        tmp[n + 1] = 1
        for i in (n + 2):(2n)
            tmp[i] = coeffs[2n - i + 2]
        end
    end

    # Apply weights and FFT
    @inbounds begin
        # Apply weights
        for i in eachindex(tmp)
            tmp[i] *= w[i]
        end

        # FFT
        fft!(tmp)
    end

    # Extract real values
    @inbounds for i in 1:n
        vals[i] = real(tmp[n - i + 1])
    end

    # Enforce symmetry if needed
    if isEven || isOdd
        half = one(TR) / 2
        @inbounds for i in 1:div(n, 2)
            j = n - i + 1
            if isEven
                s = vals[i] + vals[j]
                vals[i] = half * s
                vals[j] = half * s
            else
                d = vals[i] - vals[j]
                vals[i] = half * d
                vals[j] = -half * d
            end
        end
    end

    return vals
end

export cheb1_coeffs2vals, cheb1_coeffs2vals!, Cheb1Coeffs2ValsCache

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
