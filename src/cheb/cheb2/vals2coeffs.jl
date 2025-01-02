"""
    cheb2_vals2coeffs(vals::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}
    op::Cheb2Vals2CoeffsOp(vals::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}

Convert values at Chebyshev points of the 2nd kind into Chebyshev coefficients.

# Performance Guide
For best performance, especially in loops or repeated calls:
```julia
op = Cheb2Vals2CoeffsOp(Float64, n)
values = op(coeffs)
```

# References
- [chebfun/@chebtech2/vals2coeffs.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/vals2coeffs.m)
"""
struct Cheb2Vals2CoeffsOp{TR<:AbstractFloat,TP<:Plan}
    tmp::Vector{Complex{TR}}
    coeffs::Vector{TR}
    ifft_plan::TP

    function Cheb2Vals2CoeffsOp(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
        tmp = zeros(Complex{TR}, 2n - 2)
        coeffs = zeros(TR, n)
        ifft_plan = plan_ifft_measure!(tmp)
        return new{TR,typeof(ifft_plan)}(tmp, coeffs, ifft_plan)
    end

    function Cheb2Vals2CoeffsOp(n::TI) where {TI<:Integer}
        Cheb2Vals2CoeffsOp(Float64, n)
    end
end

function (op::Cheb2Vals2CoeffsOp{TR,TP})(
    vals::AbstractVector{TR}
) where {TR<:AbstractFloat,TP<:Plan}
    n = length(vals)

    # Trivial case
    if n <= 1
        op.coeffs .= vals
        return op.coeffs
    end

    # Determine if vals are even or odd symmetric
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

    # Mirror the values
    @inbounds for i in 1:(n - 1)
        op.tmp[i] = vals[n - i + 1]  # descending part
        op.tmp[n - 1 + i] = vals[i]  # ascending part
    end

    # Perform inverse FFT on the mirrored data
    op.ifft_plan * op.tmp

    @inbounds begin
        op.coeffs[1] = real(op.tmp[1])
        for i in 2:(n - 1)
            op.coeffs[i] = 2 * real(op.tmp[i])
        end
        op.coeffs[n] = real(op.tmp[n])
    end

    # Enforce exact symmetries
    if is_even
        @inbounds for i in 2:2:n
            op.coeffs[i] = 0
        end
    elseif is_odd
        @inbounds for i in 1:2:n
            op.coeffs[i] = 0
        end
    end

    return op.coeffs
end

function cheb2_vals2coeffs(vals::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}
    n = length(vals)
    if n <= 1
        return deepcopy(vals)
    end
    op = Cheb2Vals2CoeffsOp(TR, n)
    return op(vals)
end

export cheb2_vals2coeffs, Cheb2Vals2CoeffsOp

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

    @testset "Operator style" begin
        n = 100
        vals = rand(n)
        op = Cheb2Vals2CoeffsOp(Float64, n)

        # Test operator call
        coeffs1 = op(vals)
        coeffs2 = cheb2_vals2coeffs(vals)
        @test maximum(abs.(coeffs1 .- coeffs2)) < tol

        # Test multiple calls
        for _ in 1:10
            vals = rand(n)
            coeffs1 = op(vals)
            coeffs2 = cheb2_vals2coeffs(vals)
            @test maximum(abs.(coeffs1 .- coeffs2)) < tol
        end
    end
end