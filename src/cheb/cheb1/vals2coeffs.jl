"""
    cheb1_vals2coeffs(vals::AbstractVector{TR}) where {TR<:AbstractFloat}
    Cheb1Vals2CoeffsOp{[TR=Float64]}(n::Integer)(vals::VT) where {TR<:AbstractFloat}

Convert values at Chebyshev points of the 1st kind into Chebyshev coefficients.

# Performance Guide
For best performance, especially in loops or repeated calls:
```julia
op = Cheb1Vals2CoeffsOp{Float64}(n)
values = op(coeffs)
```

# References
- [chebfun/@chebtech1/vals2coeffs.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech1/vals2coeffs.m)
"""
struct Cheb1Vals2CoeffsOp{TR<:AbstractFloat}
    w::Vector{Complex{TR}}
    tmp::Vector{Complex{TR}}
    coeffs::Vector{Complex{TR}}
    real_coeffs::Vector{TR}
    ifft_plan::Plan{Complex{TR}}

    function Cheb1Vals2CoeffsOp{TR}(n::Integer) where {TR<:AbstractFloat}
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
        coeffs = Vector{Complex{TR}}(undef, n)
        real_coeffs = Vector{TR}(undef, n)

        # Create an inverse FFT plan with MEASURE flag for better performance
        ifft_plan = plan_ifft_measure!(tmp)

        return new{TR}(w, tmp, coeffs, real_coeffs, ifft_plan)
    end

    function Cheb1Vals2CoeffsOp(n::Integer)
        return Cheb1Vals2CoeffsOp{Float64}(n)
    end
end

function (op::Cheb1Vals2CoeffsOp)(
    vals::AbstractVector{TR}
) where {TR<:AbstractFloatOrComplex}
    type_is_float = typeisfloat(TR)

    n = length(vals)
    if n <= 1
        if type_is_float
            op.real_coeffs .= vals
            return op.real_coeffs
        else
            op.coeffs .= vals
            return op.coeffs
        end
    end

    # Check for symmetry with tolerance
    atol = type_is_float ? 10 * eps(TR) : 10 * eps(real(TR))
    isEven = true
    isOdd = true
    @inbounds for i in 1:(n ÷ 2)
        diff = abs(vals[i] - vals[n - i + 1])
        sum = abs(vals[i] + vals[n - i + 1])
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

    # Build tmp as [reverse(vals); vals] more efficiently
    @inbounds begin
        for i in 1:n
            op.tmp[i] = vals[n - i + 1]
            op.tmp[n + i] = vals[i]
        end
    end

    # Apply IFFT
    op.ifft_plan * op.tmp

    # Extract and scale coefficients
    @inbounds begin
        for k in 1:n
            op.coeffs[k] = op.tmp[k] * op.w[k]
        end
    end

    # Enforce symmetry if detected
    if isEven || isOdd
        @inbounds begin
            k_start = isEven ? 2 : 1
            op.coeffs[k_start:2:n] .= 0
        end
    end

    if type_is_float
        @inbounds for k in 1:n
            op.real_coeffs[k] = real(op.coeffs[k])
        end
        return op.real_coeffs
    else
        return op.coeffs
    end
end

function cheb1_vals2coeffs(vals::AbstractVector{TR}) where {TR<:AbstractFloatOrComplex}
    n = length(vals)
    if n <= 1
        return deepcopy(vals)
    end

    op = Cheb1Vals2CoeffsOp{real(TR)}(n)
    return op(vals)
end

export cheb1_vals2coeffs, Cheb1Vals2CoeffsOp
