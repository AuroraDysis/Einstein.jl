"""
    cheb1_vals2coeffs(values::AbstractVector{<:AbstractFloat})

Convert values at Chebyshev points of the first kind to Chebyshev coefficients.

Given a vector `values` of length `n` containing function values at Chebyshev points
of the first kind, this function returns a vector `coeffs` of Chebyshev coefficients 
such that the Chebyshev series (of the first kind) interpolates the data.

This function creates and caches computation resources internally and calls
`cheb1_vals2coeffs!` for the actual computation.

# Examples

```julia
coeffs = cheb1_vals2coeffs(values)
```

# References

- Section 4.7 of *"Chebyshev Polynomials"* by Mason & Handscomb,
  Chapman & Hall/CRC (2003).

"""
function cheb1_vals2coeffs(values::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}
    n = length(values)
    if n <= 1
        return deepcopy(values)
    end

    cache = Cheb1Vals2CoeffsCache{TR}(n)
    cheb1_vals2coeffs!(values, cache)

    return cache.coeffs
end

"""
    struct Cheb1Vals2CoeffsCache{TR}

Cache structure for `cheb1_vals2coeffs!` function to store precomputed weights,
temporary arrays, and FFT plan to speed up multiple transformations of the same
size.

Fields:

- `w::Vector{Complex{TR}}`: Weight vector of length `n`
- `tmp::Vector{Complex{TR}}`: Temporary storage vector of length `2n`
- `coeffs::Vector{TR}`: Result vector of Chebyshev coefficients of length `n`
- `plan::FFTW.cFFTWPlan`: FFTW plan for inverse FFT of size `2n`

"""
struct Cheb1Vals2CoeffsCache{TR<:AbstractFloat,TP<:Plan}
    w::Vector{Complex{TR}}
    tmp::Vector{Complex{TR}}
    coeffs::Vector{TR}
    ifft_plan::TP

    function Cheb1Vals2CoeffsCache{TR}(n::Integer) where {TR<:AbstractFloat}
        # Precompute weights
        w = Vector{Complex{TR}}(undef, n)
        @inbounds begin
            im_pi_over_2n = im * convert(TR, π) / (2n)
            for k in 0:(n - 1)
                w[k + 1] = 2 * exp(k * im_pi_over_2n)
            end
            w[1] /= 2  # Special case for k=0
        end

        # Prepare temporary array for FFT
        tmp = Vector{Complex{TR}}(undef, 2n)

        coeffs = Vector{TR}(undef, n)

        # Create an inverse FFT plan with MEASURE flag for better performance
        ifft_plan = plan_ifft_measure!(tmp)

        return new{TR,typeof(ifft_plan)}(w, tmp, coeffs, ifft_plan)
    end
end

"""
    cheb1_vals2coeffs!(values::AbstractVector{<:AbstractFloat},
                       cache::Cheb1Vals2CoeffsCache{TR}) where {TR}

In-place version of `cheb1_vals2coeffs` that uses a cache for efficiency.

Computes the Chebyshev coefficients corresponding to the given `values` at
Chebyshev points of the first kind, storing the result in `cache.coeffs`.

# Arguments

- `values`: Vector of function values at Chebyshev points of the first kind.
- `cache`: An instance of `Cheb1Vals2CoeffsCache` for storing precomputed data
           and temporary arrays.

# Returns

- The vector `cache.coeffs`, containing the Chebyshev coefficients.

"""
function cheb1_vals2coeffs!(
    values::VT, cache::Cheb1Vals2CoeffsCache{TR}
) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}
    n = length(values)
    if n <= 1
        cache.coeffs .= values
        return cache.coeffs
    end

    # Check for symmetry with tolerance
    atol = 10 * eps(TR)
    isEven = true
    isOdd = true
    @inbounds for i in 1:(n ÷ 2)
        diff = abs(values[i] - values[n - i + 1])
        sum = abs(values[i] + values[n - i + 1])
        if diff > atol
            isEven = false
        end
        if sum > atol
            isOdd = false
        end
        if !isEven && !isOdd
            break
        end
    end

    tmp = cache.tmp
    w = cache.w
    coeffs = cache.coeffs
    ifft_plan = cache.ifft_plan

    # Build tmp as [reverse(values); values] more efficiently
    @inbounds begin
        for i in 1:n
            tmp[i] = Complex{TR}(values[n - i + 1])
            tmp[n + i] = Complex{TR}(values[i])
        end
    end

    # Apply IFFT
    ifft_plan * tmp

    # Extract and scale coefficients
    @inbounds begin
        for k in 1:n
            coeffs[k] = real(tmp[k] * w[k])
        end
    end

    # Enforce symmetry if detected
    if isEven || isOdd
        @inbounds begin
            k_start = isEven ? 2 : 1
            coeffs[k_start:2:n] .= 0
        end
    end

    return coeffs
end

export cheb1_vals2coeffs, cheb1_vals2coeffs!, Cheb1Vals2CoeffsCache

@testset "Cheb1 vals2coeffs tests" begin
    tol = 100 * eps()

    @testset "Single value conversion" begin
        v = sqrt(2)
        c = cheb1_vals2coeffs([v])
        @test v ≈ c[1]
    end

    @testset "Even case tests" begin
        v = Float64[1:6;]
        cTrue = [
            7 / 2,
            sqrt(6) / 2 + 5 * sqrt(2) / 6,
            0,
            sqrt(2) / 6,
            0,
            sqrt(6) / 2 - 5 * sqrt(2) / 6,
        ]
        c = cheb1_vals2coeffs(v)
        @test norm(c - cTrue, Inf) < tol
        @test all(x -> abs(imag(x)) < tol, c)
    end

    @testset "Odd case tests" begin
        v = Float64[1:5;]
        cTrue = [
            3,
            (2 / 5) * (sqrt((5 - sqrt(5)) / 2) + 2 * sqrt((5 + sqrt(5)) / 2)),
            0,
            (2 / 5) * (2 * sqrt((5 - sqrt(5)) / 2) - sqrt((5 + sqrt(5)) / 2)),
            0,
        ]
        c = cheb1_vals2coeffs(v)
        @test norm(c - cTrue, Inf) < tol
        @test all(x -> abs(imag(x)) < tol, c)
    end

    @testset "Symmetry preservation" begin
        v = kron([1 -1; 1 1], ones(10, 1))
        c1 = cheb1_vals2coeffs(@view(v[:, 1]))
        c2 = cheb1_vals2coeffs(@view(v[:, 2]))
        @test norm(c1[2:2:end], Inf) ≈ 0
        @test norm(c2[1:2:end], Inf) ≈ 0
    end
end