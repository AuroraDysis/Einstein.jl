"""
    cheb2_coeffs2vals(coeffs::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}

Convert Chebyshev coefficients of the 2nd kind to values at Chebyshev points.

# Arguments
- `coeffs::AbstractVector{<:AbstractFloat}`: Vector of Chebyshev coefficients

# Returns
- Vector of values at Chebyshev points of the 2nd kind

# Description
This function transforms Chebyshev coefficients to values at Chebyshev points using an FFT-based
algorithm. For a polynomial of degree n-1, the transformation maps n coefficients to n values
at the Chebyshev points of the 2nd kind: cos(jπ/(n-1)) for j = 0:n-1.

# Example
```julia
coeffs = [1.0, 2.0, 3.0]  # Chebyshev coefficients
vals = cheb2_coeffs2vals(coeffs)  # Values at Chebyshev points
```
"""
function cheb2_coeffs2vals(coeffs::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}
    n = length(coeffs)

    # Trivial case (constant or empty)
    if n <= 1
        return deepcopy(coeffs)
    end

    op = Cheb2Coeffs2ValsOp(TR, n)
    return op(coeffs)
end

"""
    Cheb2Coeffs2ValsOp{TR<:AbstractFloat,TP<:Plan}

A pre-planned operator for efficient conversion from Chebyshev coefficients to values.

# Fields
- `tmp::Vector{Complex{TR}}`: Temporary storage for FFT computation
- `vals::Vector{TR}`: Storage for the resulting values
- `fft_plan::TP`: Pre-computed FFT plan for efficient transformation

This struct provides a reusable transformation operator that avoids repeated memory allocation
and FFT plan creation when performing multiple transformations of the same size.
"""
struct Cheb2Coeffs2ValsOp{TR<:AbstractFloat,TP<:Plan}
    tmp::Vector{Complex{TR}}
    vals::Vector{TR}
    fft_plan::TP

    function Cheb2Coeffs2ValsOp(::Type{TR}, n::TI) where {TI<:Integer,TR<:AbstractFloat}
        tmp = zeros(Complex{TR}, 2n - 2)
        vals = zeros(TR, n)
        fft_plan = plan_fft_measure!(tmp)
        return new{TR,typeof(fft_plan)}(tmp, vals, fft_plan)
    end

    function Cheb2Coeffs2ValsOp(n::TI) where {TI<:Integer}
        return Cheb2Coeffs2ValsOp(Float64, n)
    end
end

"""
    (op::Cheb2Coeffs2ValsOp{TR,TP})(coeffs::VT)

Apply the Chebyshev coefficient to values transformation using a pre-planned operator.

# Arguments
- `coeffs::AbstractVector{<:AbstractFloat}`: Vector of Chebyshev coefficients

# Returns
- Vector of values at Chebyshev points

# Implementation Details
1. For n coefficients, creates a padded array of length 2n-2
2. Applies symmetry to create the correct even/odd extension
3. Performs FFT transformation
4. Extracts and scales the real parts to obtain final values
5. Enforces symmetry properties based on coefficient pattern:
   - Even symmetry if even-indexed coefficients are zero
   - Odd symmetry if odd-indexed coefficients are zero
"""
function (op::Cheb2Coeffs2ValsOp{TR,TP})(coeffs::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR},TP<:Plan}
    n = length(coeffs)
    if n <= 1
        op.vals .= coeffs
        return op.vals
    end

    vals = op.vals
    tmp = op.tmp
    fft_plan = op.fft_plan

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
        fft_plan * tmp

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

export cheb2_coeffs2vals, Cheb2Coeffs2ValsOp

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
