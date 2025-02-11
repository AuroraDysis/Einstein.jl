"""
    cheb2_vals2coeffs(vals::AbstractVector{TF}) where {TF<:AbstractFloat}
    ChebyshevSecondKindAnalysis{[TF=Float64]}(n::Integer)(vals::AbstractVector{TF}) where {TF<:AbstractFloat}

Convert values at Chebyshev points of the 2nd kind into Chebyshev coefficients.

# Performance Guide
For best performance, especially in loops or repeated calls:
```julia
op = ChebyshevSecondKindAnalysis{Float64}(n)
values = op(coeffs)
```

# References
- [chebfun/@chebtech2/vals2coeffs.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/vals2coeffs.m)
"""
struct ChebyshevSecondKindAnalysis{TF<:AbstractFloat} <:
       AbstractChebyshevAnalysisImplementation
    tmp::Vector{Complex{TF}}
    coeffs::Vector{Complex{TF}}
    real_coeffs::Vector{TF}
    ifft_plan::Plan{Complex{TF}}

    function ChebyshevSecondKindAnalysis{TF}(n::Integer) where {TF<:AbstractFloat}
        tmp = zeros(Complex{TF}, 2n - 2)
        coeffs = zeros(Complex{TF}, n)
        real_coeffs = zeros(TF, n)
        ifft_plan = plan_ifft_measure!(tmp)
        return new{TF}(tmp, coeffs, real_coeffs, ifft_plan)
    end

    function ChebyshevSecondKindAnalysis(n::Integer)
        return ChebyshevSecondKindAnalysis{Float64}(n)
    end
end

function (op::ChebyshevSecondKindAnalysis{TF})(
    vals::AbstractVector{TFC}
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}}
    type_is_float = typeisfloat(TFC)

    n = length(vals)

    # Trivial case
    if n <= 1
        if type_is_float
            op.real_coeffs .= vals
            return op.real_coeffs
        else
            op.coeffs .= vals
            return op.coeffs
        end
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

    # Scale the interior coefficients
    @inbounds begin
        op.coeffs[1] = op.tmp[1]
        for i in 2:(n - 1)
            op.coeffs[i] = 2 * op.tmp[i]
        end
        op.coeffs[n] = op.tmp[n]
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

    if type_is_float
        @inbounds for k in 1:n
            op.real_coeffs[k] = real(op.coeffs[k])
        end
        return op.real_coeffs
    else
        return op.coeffs
    end
end

function cheb_analysis(
    ::ChebyshevSecondKindNode, ::Type{TF}, n::Integer
) where {TF<:AbstractFloat}
    return ChebyshevSecondKindAnalysis{TF}(n)
end

function cheb2_vals2coeffs(vals::AbstractVector{TFC}) where {TFC<:AbstractFloatOrComplex}
    n = length(vals)
    if n <= 1
        return deepcopy(vals)
    end
    op = ChebyshevSecondKindAnalysis{real(TFC)}(n)
    return op(vals)
end

export ChebyshevSecondKindAnalysis, cheb2_vals2coeffs
