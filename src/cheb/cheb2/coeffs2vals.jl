using FFTW: fft, ifft
using FastTransforms: fft, ifft

"""
    cheb2_coeffs2vals(coeffs::AbstractVector{<:AbstractFloat}) -> Vector

Convert Chebyshev coefficients of the second kind to values at Chebyshev nodes using a DCT-I equivalent
transform.

# Arguments
- `coeffs::AbstractVector{<:AbstractFloat}`: Vector of Chebyshev coefficients in descending order

# Returns
- Vector of values at Chebyshev nodes of the second kind: cos(jπ/n), j = 0,1,...,n

# Mathematical Background
The function implements the transform from coefficient space to physical space for Chebyshev
polynomials of the second kind Uₙ(x). The transformation preserves symmetries:
- Even coefficients map to even functions
- Odd coefficients map to odd functions

# Implementation Details
1. Scales interior coefficients by 1/2
2. Mirrors coefficients and applies FFT
3. Enforces even/odd symmetries in the result

# Examples
```julia
# Convert a constant function
julia> cheb2_coeffs2vals([1.0])
1-element Vector{Float64}:
 1.0

# Convert linear coefficients
julia> cheb2_coeffs2vals([1.0, 2.0])
2-element Vector{Float64}:
  3.0
 -1.0
```

# References
1. Trefethen, L. N. (2000). Spectral Methods in MATLAB. SIAM.
2. Boyd, J. P. (2001). Chebyshev and Fourier Spectral Methods. Dover.

See also: [`cheb1_coeffs2vals`](@ref)
"""
function cheb2_coeffs2vals(coeffs::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}
    n = length(coeffs)

    # Trivial case (constant or empty)
    if n <= 1
        return deepcopy(coeffs)
    end

    # Determine which columns are purely even or purely odd based on middle coefficients
    isEven = all(x -> x ≈ 0, @view(coeffs[2:2:end]))
    isOdd = all(x -> x ≈ 0, @view(coeffs[1:2:end]))

    # Scale the interior rows by 1/2
    half = one(TR) / 2
    coeffs[2:(end - 1), :] .*= half

    # Mirror the coefficients for a DCT-I using FFT
    tmp = vcat(coeffs, coeffs[(n - 1):-1:2, :])

    # apply the FFT.
    values = real.(fft(tmp))

    # Flip and truncate to size n
    values = values[n:-1:1]

    # Enforce symmetry in each column based on isEven and isOdd flags
    if isEven
        # (values + flipped(values)) / 2
        values = half .* (values .+ reverse(values))
    elseif isOdd
        # (values - flipped(values)) / 2
        values = half .* (values .- reverse(values))
    end

    return values
end

export cheb2_coeffs2vals

@testset "cheb2_coeffs2vals" begin
    using LinearAlgebra

    # Set tolerance
    tol = 100 * eps()

    # Test single coefficient conversion
    c = [sqrt(2)]
    v = cheb2_coeffs2vals(c)
    @test c ≈ v

    # Test simple data conversion
    c = collect(1.0:5.0)
    # Exact values
    vTrue = [3; -4 + sqrt(2); 3; -4 - sqrt(2); 15]

    # Test real branch
    v = cheb2_coeffs2vals(c)
    @test norm(v - vTrue, Inf) < tol
    @test all(iszero, imag.(v))

    # Test symmetry preservation
    c = kron(ones(10), Matrix{Float64}(I, 2, 2))
    c1 = @view c[:, 1]
    c2 = @view c[:, 2]
    v1 = cheb2_coeffs2vals(c1)
    v2 = cheb2_coeffs2vals(c2)
    @test norm(v1 - reverse(v1), Inf) ≈ 0
    @test norm(v2 + reverse(v2), Inf) ≈ 0
end
