"""
    cheb1_vals2coeffs(vals::AbstractVector{TF}) where {TF<:AbstractFloat}
    ChebyshevFirstKindAnalysis{[TF=Float64]}(n::Integer)(vals::VT) where {TF<:AbstractFloat}

Convert values at Chebyshev points of the 1st kind into Chebyshev coefficients.

# Performance Guide
For best performance, especially in loops or repeated calls:
```julia
op = ChebyshevFirstKindAnalysis{Float64}(n)
values = op(coeffs)
```

# References
- [chebfun/@chebtech1/vals2coeffs.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech1/vals2coeffs.m)
"""
struct ChebyshevFirstKindAnalysis{TF<:AbstractFloat} <:
       AbstractChebyshevAnalysisImplementation
    w::Vector{Complex{TF}}
    tmp::Vector{Complex{TF}}
    coeffs::Vector{Complex{TF}}
    real_coeffs::Vector{TF}
    ifft_plan::Plan{Complex{TF}}

    function ChebyshevFirstKindAnalysis{TF}(n::Integer) where {TF<:AbstractFloat}
        # Precompute weights
        w = Vector{Complex{TF}}(undef, n)
        @inbounds begin
            im_pi_over_2n = im * convert(TF, π) / (2n)
            for k in 0:(n - 1)
                w[k + 1] = 2 * exp(k * im_pi_over_2n)
            end
            w[1] /= 2  # Special case for k=0
        end

        # Prepare temporary array for FFT
        tmp = Vector{Complex{TF}}(undef, 2n)
        coeffs = Vector{Complex{TF}}(undef, n)
        real_coeffs = Vector{TF}(undef, n)

        # Create an inverse FFT plan with MEASURE flag for better performance
        ifft_plan = plan_ifft_measure!(tmp)

        return new{TF}(w, tmp, coeffs, real_coeffs, ifft_plan)
    end

    function ChebyshevFirstKindAnalysis(n::Integer)
        return ChebyshevFirstKindAnalysis{Float64}(n)
    end
end

function (op::ChebyshevFirstKindAnalysis{TF})(
    vals::AbstractVector{TFC}
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}}
    type_is_float = typeisfloat(TFC)

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
    atol = type_is_float ? 10 * eps(TF) : 10 * eps(real(TF))
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

function _cheb_analysis(
    ::ChebyshevT, ::Type{TF}, n::Integer
) where {TF<:AbstractFloat}
    return ChebyshevFirstKindAnalysis{TF}(n)
end

function cheb1_vals2coeffs(vals::AbstractVector{TFC}) where {TFC<:AbstractFloatOrComplex}
    n = length(vals)
    if n <= 1
        return deepcopy(vals)
    end

    op = ChebyshevFirstKindAnalysis{real(TFC)}(n)
    return op(vals)
end

export ChebyshevFirstKindAnalysis, cheb1_vals2coeffs
