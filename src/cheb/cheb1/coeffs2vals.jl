"""
    cheb1_coeffs2vals(coeffs::AbstractVector{TR}) where {TR<:AbstractFloat}
    Cheb1Coeffs2ValsOp{[TR=Float64]}(n::Integer)(coeffs::AbstractVector{TR}) where {TR<:AbstractFloat}

Convert Chebyshev coefficients to values at Chebyshev points of the 1st kind.

# Performance Guide
For best performance, especially in loops or repeated calls:
```julia
op = Cheb1Coeffs2ValsOp{Float64}(n)
values = op(coeffs)
```

# References
- [chebfun/@chebtech1/coeffs2vals.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech1/coeffs2vals.m)
"""
struct Cheb1Coeffs2ValsOp{TR<:AbstractFloat}
    w::Vector{Complex{TR}}    # Weight vector
    tmp::Vector{Complex{TR}}  # Temporary storage
    vals::Vector{TR}          # Result storage
    fft_plan::Plan{Complex{TR}}        # fft plan

    function Cheb1Coeffs2ValsOp{TR}(n::Integer) where {TR<:AbstractFloat}
        # Precompute weights
        w = Vector{Complex{TR}}(undef, 2n)
        @inbounds begin
            half = one(TR) / 2
            m_im_pi_over_2n = -im * convert(TR, π) / (2n)
            for k in 0:(n - 1)
                w[k + 1] = exp(k * m_im_pi_over_2n) * half
            end
            w[1] *= 2
            w[n + 1] = 0
            for k in (n + 1):(2n - 1)
                w[k + 1] = -exp(k * m_im_pi_over_2n) * half
            end
        end
        tmp = Vector{Complex{TR}}(undef, 2n)
        vals = Vector{TR}(undef, n)
        fft_plan = plan_fft_measure!(tmp)
        return new{TR}(w, tmp, vals, fft_plan)
    end

    function Cheb1Coeffs2ValsOp(n::Integer)
        return Cheb1Coeffs2ValsOp{Float64}(n)
    end
end

function (op::Cheb1Coeffs2ValsOp{TR})(coeffs::AbstractVector{TR}) where {TR<:AbstractFloat}
    n = length(coeffs)
    if n <= 1
        op.vals .= coeffs
        return op.vals
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

    # Extract real values
    @inbounds for i in 1:n
        vals[i] = real(tmp[n - i + 1])
    end

    # Enforce symmetry if needed
    if isEven || isOdd
        half = one(TR) / 2
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

    return vals
end

function cheb1_coeffs2vals(coeffs::AbstractVector{TR}) where {TR<:AbstractFloat}
    n = length(coeffs)
    if n <= 1
        return deepcopy(coeffs)
    end

    op = Cheb1Coeffs2ValsOp{TR}(n)
    return op(coeffs)
end

export cheb1_coeffs2vals, Cheb1Coeffs2ValsOp
