"""
    cheb2_coeffs2vals(coeffs::AbstractVector{<:AbstractFloat}) -> Vector
    cheb2_coeffs2vals!(coeffs, cache::Cheb2Coeffs2ValsCache) -> Vector

Convert Chebyshev coefficients of the second kind to values at Chebyshev nodes using a DCT-I equivalent
transform.

# Performance Guide
For best performance, especially in loops or repeated calls:
1. Create a cache: `cache = Cheb2Coeffs2ValsCache{Float64}(n)`
2. Use the in-place version: `cheb2_coeffs2vals!(coeffs, cache)`

This avoids repeated memory allocations and can be significantly faster.

# Examples
```julia
# Single transform (allocating version)
julia> v = cheb2_coeffs2vals([1.0, 2.0])

# Multiple transforms (cached version for better performance)
julia> cache = Cheb2Coeffs2ValsCache{Float64}(2)
julia> for coeffs in coefficient_arrays
           v = cheb2_coeffs2vals!(coeffs, cache)
           # ... use v ...
       end
```

# Arguments
- `coeffs::AbstractVector{<:AbstractFloat}`: Vector of Chebyshev coefficients in descending order
- `cache::Cheb2Coeffs2ValsCache`: Pre-allocated workspace (for in-place version)

# Returns
- Vector of values at Chebyshev nodes of the second kind: [`cheb2_pts`](@ref)

# Cache Creation
```julia
# Create a cache for size n transforms
cache = Cheb2Coeffs2ValsCache{Float64}(n)
```

# Mathematical Background
The function implements the transform from coefficient space to physical space for Chebyshev
polynomials of the second kind \$U_n(x)\$. The transformation preserves symmetries:
- Even coefficients map to even functions
- Odd coefficients map to odd functions

# Implementation Details
1. Scales interior coefficients by 1/2
2. Mirrors coefficients and applies FFT
3. Enforces even/odd symmetries in the result

# References
1. Trefethen, L. N. (2000). Spectral Methods in MATLAB. SIAM.
2. Boyd, J. P. (2001). Chebyshev and Fourier Spectral Methods. Dover.

See also: [`cheb1_coeffs2vals`](@ref), [`Cheb2Coeffs2ValsCache`](@ref)
"""
function cheb2_coeffs2vals(coeffs::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}
    n = length(coeffs)

    # Trivial case (constant or empty)
    if n <= 1
        return deepcopy(coeffs)
    end

    cache = Cheb2Coeffs2ValsCache{TR}(n)
    cheb2_coeffs2vals!(coeffs, cache)

    return cache.vals
end

"""
    Cheb2Coeffs2ValsCache{T}

Pre-allocated workspace for Chebyshev coefficient to values transformation.
Using this cache can significantly improve performance when performing multiple transforms
of the same size.

# Fields
- `tmp::Vector{Complex{T}}`: Temporary storage for FFT computation
- `vals::Vector{T}`: Storage for the final result

# Example
```julia
# Create cache for size 100 transforms
cache = Cheb2Coeffs2ValsCache{Float64}(100)

# Use cache repeatedly
for i in 1:1000
    vals = cheb2_coeffs2vals!(coeffs, cache)
end
```
"""
struct Cheb2Coeffs2ValsCache{TR}
    tmp::Vector{Complex{TR}}
    vals::Vector{TR}

    function Cheb2Coeffs2ValsCache{TR}(n::TI) where {TI<:Integer,TR<:AbstractFloat}
        return new(zeros(Complex{TR}, 2n - 2), zeros(TR, n))
    end
end

function cheb2_coeffs2vals!(
    coeffs::VT, cache::Cheb2Coeffs2ValsCache{TR}
) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}
    n = length(coeffs)
    if n <= 1
        cache.vals .= coeffs
        return cache.vals
    end

    vals = cache.vals
    tmp = cache.tmp

    # Determine which columns are purely even or purely odd based on middle coefficients
    isEven = all(x -> x ≈ 0, @view(coeffs[2:2:end]))
    isOdd = all(x -> x ≈ 0, @view(coeffs[1:2:end]))

    half = one(TR) / 2
    @inbounds begin
        tmp[1] = coeffs[1]
        for i in 2:(n - 1)
            hc = half * coeffs[i]
            tmp[i] = hc
            tmp[2n - i] = hc
        end
        tmp[n] = coeffs[n]

        # FFT into vals
        fft!(tmp)

        # Flip/truncate inside vals
        for i in 1:n
            vals[i] = real(tmp[n - i + 1])
        end
    end

    # In-place symmetry enforcement (reuse logic from original):
    if isEven
        @inbounds for i in 1:div(length(vals), 2)
            j = length(vals) - i + 1
            s = vals[i] + vals[j]
            vals[i] = half * s
            vals[j] = half * s
        end
    elseif isOdd
        @inbounds for i in 1:div(length(vals), 2)
            j = length(vals) - i + 1
            d = vals[i] - vals[j]
            vals[i] = half * d
            vals[j] = -half * d
        end
    end

    return vals
end

export cheb2_coeffs2vals, cheb2_coeffs2vals!, Cheb2Coeffs2ValsCache

@testset "cheb2_coeffs2vals" begin
    using LinearAlgebra

    # Set tolerance
    tol = 100 * eps()

    # Test single coefficient conversion
    c = [sqrt(2)]
    v = cheb2_coeffs2vals(c)
    @test c ≈ v

    c = collect(1.0:5.0)
    vTrue = [3; -4 + sqrt(2); 3; -4 - sqrt(2); 15]
    v = cheb2_coeffs2vals(c)
    @test norm(v - vTrue, Inf) < tol

    c = collect(1.0:6.0)
    vTrue = [-3; 7 / 2; -(11 / 2) + sqrt(5); 7 / 2; -(11 / 2) - sqrt(5); 21]
    v = cheb2_coeffs2vals(c)
    @test norm(v - vTrue, Inf) < tol

    # Test symmetry preservation
    c = kron(ones(10), Matrix{Float64}(I, 2, 2))
    c1 = @view c[:, 1]
    c2 = @view c[:, 2]
    v1 = cheb2_coeffs2vals(c1)
    v2 = cheb2_coeffs2vals(c2)
    @test norm(v1 - reverse(v1), Inf) ≈ 0
    @test norm(v2 + reverse(v2), Inf) ≈ 0
end
