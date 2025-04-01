"""
    chebtech1_coeffs2vals(coeffs::AbstractVector{TFC}) where {TFC<:AbstractFloatOrComplex}
    ChebyshevFirstKindSynthesis{[TF=Float64]}(n::Integer)(coeffs::AbstractVector{TFC})

Convert Chebyshev coefficients to values at Chebyshev points of the 1st kind.

# Performance Guide
For best performance, especially in loops or repeated calls:
```julia
op = ChebyshevFirstKindSynthesis{Float64}(n)
values = op(coeffs)
```

# References
- [chebfun/@chebtech1/coeffs2vals.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech1/coeffs2vals.m)
"""
struct ChebyshevFirstKindSynthesis{TF<:AbstractFloat} <:
       AbstractChebyshevSynthesisImplementation
    w::Vector{Complex{TF}}    # Weight vector
    tmp::Vector{Complex{TF}}  # Temporary storage
    vals::Vector{Complex{TF}} # values
    real_vals::Vector{TF} # values
    fft_plan::Plan{Complex{TF}}        # fft plan

    function ChebyshevFirstKindSynthesis{TF}(n::Integer) where {TF<:AbstractFloat}
        # Precompute weights
        w = Vector{Complex{TF}}(undef, 2n)
        @inbounds begin
            half = one(TF) / 2
            m_im_pi_over_2n = -im * convert(TF, π) / (2n)
            for k in 0:(n - 1)
                w[k + 1] = exp(k * m_im_pi_over_2n) * half
            end
            w[1] *= 2
            w[n + 1] = 0
            for k in (n + 1):(2n - 1)
                w[k + 1] = -exp(k * m_im_pi_over_2n) * half
            end
        end
        tmp = Vector{Complex{TF}}(undef, 2n)
        vals = Vector{Complex{TF}}(undef, n)
        real_vals = Vector{TF}(undef, n)
        fft_plan = plan_fft_measure!(tmp)
        return new{TF}(w, tmp, vals, real_vals, fft_plan)
    end

    function ChebyshevFirstKindSynthesis(n::Integer)
        return ChebyshevFirstKindSynthesis{Float64}(n)
    end
end

function (op::ChebyshevFirstKindSynthesis{TF})(
    coeffs::AbstractVector{TFC}
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}}
    type_is_float = typeisfloat(TFC)

    n = length(coeffs)
    if n <= 1
        if type_is_float
            op.real_vals .= coeffs
            return op.real_vals
        else
            op.vals .= coeffs
            return op.vals
        end
    end

    w = op.w
    tmp = op.tmp
    vals = op.vals
    fft_plan = op.fft_plan

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
        fft_plan * tmp
    end

    # Extract values
    @inbounds for i in 1:n
        vals[i] = tmp[n - i + 1]
    end

    # Enforce symmetry if needed
    if isEven || isOdd
        half = one(TF) / 2
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

    if type_is_float
        @inbounds for k in 1:n
            op.real_vals[k] = real(op.vals[k])
        end
        return op.real_vals
    else
        return op.vals
    end
end

function _cheb_synthesis(
    ::ChebyshevT, ::Type{TF}, n::Integer
) where {TF<:AbstractFloat}
    return ChebyshevFirstKindSynthesis{TF}(n)
end

function chebtech1_coeffs2vals(coeffs::AbstractVector{TFC}) where {TFC<:AbstractFloatOrComplex}
    n = length(coeffs)

    if n <= 1
        return deepcopy(coeffs)
    end

    op = ChebyshevFirstKindSynthesis{real(TFC)}(n)
    return op(coeffs)
end

export ChebyshevFirstKindSynthesis, chebtech1_coeffs2vals
