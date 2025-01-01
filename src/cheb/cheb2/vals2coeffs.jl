"""
    cheb2_vals2coeffs(vals::AbstractVector{T}) where T<:AbstractFloat

Convert values sampled at Chebyshev points of the second kind into their corresponding
Chebyshev coefficients. This function allocates a cache internally to perform the 
inverse Discrete Cosine Transform of Type I (mirrored FFT) and returns a new array 
containing the Chebyshev coefficients.

# Description

Given an input vector `vals` of length `n` representing function values at Chebyshev points
of the second kind, `cheb2_vals2coeffs` computes the Chebyshev coefficients `c` such that:

    
    f(x) = c[1]*T₀(x) + c[2]*T₁(x) + ... + c[n]*Tₙ₋₁(x)
    

where Tₖ(x) are the Chebyshev polynomials of the first kind. Internally, this function:

1. Detects trivial cases (e.g. `n <= 1`).
2. Constructs a mirror of the input values to emulate an inverse DCT using FFT.
3. Performs an inverse FFT and rescales the interior coefficients by 2.
4. Enforces exact symmetries (even/odd) where detected.

# Example
```julia
vals = [0.0, 1.0, 2.0, 1.0, 0.0]  # Example of 5 sample values
coeffs = cheb2_vals2coeffs(vals)
```
"""
function cheb2_vals2coeffs(vals::AbstractVector{T}) where {T<:AbstractFloat}
    n = length(vals)

    # Trivial case (constant or empty)
    if n <= 1
        return deepcopy(vals)
    end

    # Create a cache for repeated transforms
    cache = Cheb2Vals2CoeffsCache{T}(n)
    return cheb2_vals2coeffs!(vals, cache)
end

"""
    Cheb2Vals2CoeffsCache{T}

Pre-allocated workspace for Chebyshev values-to-coefficients transformations.
Using this cache can significantly improve performance when performing multiple
transforms of the same size.

# Fields
- `tmp::Vector{Complex{T}}`: Temporary storage for the mirrored FFT computation
- `coeffs::Vector{T}`: Storage for the final result (the Chebyshev coefficients)

# Example
```julia
# Create cache for size-100 transforms
cache = Cheb2Vals2CoeffsCache{Float64}(100)

# Reuse the cache for multiple transforms
for i in 1:1000
    some_new_vals = rand(100)  # or your own data
    coeffs = cheb2_vals2coeffs!(some_new_vals, cache)
end
```
"""
struct Cheb2Vals2CoeffsCache{TR<:AbstractFloat,TP<:Plan}
    tmp::Vector{Complex{TR}}
    coeffs::Vector{TR}
    ifft_plan::TP

    function Cheb2Vals2CoeffsCache{TR}(n::Integer) where {TR<:AbstractFloat}
        tmp = zeros(Complex{TR}, 2n - 2)
        coeffs = zeros(TR, n)
        ifft_plan = plan_ifft_measure!(tmp)
        return new{TR,typeof(ifft_plan)}(tmp, coeffs, ifft_plan)
    end
end

function cheb2_vals2coeffs!(
    vals::AbstractVector{T}, cache::Cheb2Vals2CoeffsCache{T}
) where {T<:AbstractFloat}
    n = length(vals)

    # Trivial case
    if n <= 1
        cache.coeffs .= vals
        return cache.coeffs
    end

    # Determine if vals are even or odd symmetric: 
    # Compare vals with its reverse to see if they are negatives or equal.
    # Use explicit loop to avoid allocation and compute max difference
    is_even = true
    is_odd = true
    @inbounds for i in 1:(n ÷ 2)
        diff = abs(vals[i] - vals[n - i + 1])
        sum = abs(vals[i] + vals[n - i + 1])
        if !(diff ≈ 0)
            is_even = false
        end
        if !(sum ≈ 0)
            is_odd = false
        end
        # Early exit if neither symmetry is possible
        if !is_even && !is_odd
            break
        end
    end

    tmp = cache.tmp
    coeffs = cache.coeffs
    ifft_plan = cache.ifft_plan

    # Mirror the values (similar to MATLAB's [vals(n:-1:2) ; vals(1:n-1)])
    # The mirrored data occupies tmp[1 : 2n-2].
    @inbounds for i in 1:(n - 1)
        tmp[i] = vals[n - i + 1]  # descending part
        tmp[n - 1 + i] = vals[i]  # ascending part
    end

    # Perform inverse FFT on the mirrored data
    ifft_plan * tmp

    @inbounds begin
        coeffs[1] = real(tmp[1])
        for i in 2:(n - 1)
            coeffs[i] = 2 * real(tmp[i])
        end
        coeffs[n] = real(tmp[n])
    end

    # Enforce exact symmetries (if the original data is purely even or purely odd)
    if is_even
        # Zero out the odd coefficients
        @inbounds for i in 2:2:n
            coeffs[i] = 0
        end
    elseif is_odd
        # Zero out the even coefficients
        @inbounds for i in 1:2:n
            coeffs[i] = 0
        end
    end

    return coeffs
end

export cheb2_vals2coeffs, cheb2_vals2coeffs!, Cheb2Vals2CoeffsCache

@testset "cheb2_vals2coeffs" begin
    # Set tolerance
    tol = 100 * eps()

    @testset "Single value" begin
        v = sqrt(2)
        c = cheb2_vals2coeffs([v])
        @test v ≈ c[1]
    end

    @testset "Simple data" begin
        v = collect(1.0:5.0)
        # Exact coefficients
        cTrue = [3.0, 1 + 1 / sqrt(2), 0.0, 1 - 1 / sqrt(2), 0.0]
        c = cheb2_vals2coeffs(v)
        @test maximum(abs.(c .- cTrue)) < tol
        @test all(abs.(imag.(c)) .== 0)
    end

    @testset "Array input" begin
        v = collect(1.0:5.0)
        cTrue = [3.0, 1 + 1 / sqrt(2), 0.0, 1 - 1 / sqrt(2), 0.0]

        # Test forward and reversed arrays
        c1 = cheb2_vals2coeffs(v)
        c2 = cheb2_vals2coeffs(reverse(v))

        tmp = ones(length(cTrue))
        tmp[(end - 1):-2:1] .= -1

        @test maximum(abs.(c1 .- cTrue)) < tol
        @test maximum(abs.(c2 .- (tmp .* cTrue))) < tol
    end

    @testset "Symmetry preservation" begin
        # Create test data with even/odd symmetry
        n = 10
        v_even = repeat([1.0], n)
        v_odd = repeat([-1.0, 1.0], n ÷ 2)

        c_even = cheb2_vals2coeffs(v_even)
        c_odd = cheb2_vals2coeffs(v_odd)

        # Even coefficients should have zero odd terms
        @test all(abs.(c_even[2:2:end]) .< tol)
        # Odd coefficients should have zero even terms
        @test all(abs.(c_odd[1:2:end]) .< tol)
    end

    @testset "Cache reuse" begin
        n = 100
        cache = Cheb2Vals2CoeffsCache{Float64}(n)
        v = rand(n)

        # Results should be the same with and without cache
        c1 = cheb2_vals2coeffs(v)
        c2 = cheb2_vals2coeffs!(v, cache)

        @test c1 ≈ c2
    end
end