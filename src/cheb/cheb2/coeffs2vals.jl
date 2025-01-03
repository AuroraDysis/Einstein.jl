"""
    cheb2_coeffs2vals(coeffs::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}
    op::Cheb2Coeffs2ValsOp([TR=Float64], n::Integer)(coeffs::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}

Convert Chebyshev coefficients to values at Chebyshev points of the 2nd kind.

# Performance Guide
For best performance, especially in loops or repeated calls:
```julia
op = Cheb2Coeffs2ValsOp(Float64, n)
values = op(coeffs)
```

# References
- [chebfun/@chebtech2/coeffs2vals.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/coeffs2vals.m)
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

struct Cheb2Coeffs2ValsOp{TR<:AbstractFloat,TP<:Plan}
    tmp::Vector{Complex{TR}}
    vals::Vector{TR}
    fft_plan::TP

    function Cheb2Coeffs2ValsOp(::Type{TR}, n::Integer) where {TR<:AbstractFloat}
        tmp = zeros(Complex{TR}, 2n - 2)
        vals = zeros(TR, n)
        fft_plan = plan_fft_measure!(tmp)
        return new{TR,typeof(fft_plan)}(tmp, vals, fft_plan)
    end

    function Cheb2Coeffs2ValsOp(n::Integer)
        return Cheb2Coeffs2ValsOp(Float64, n)
    end
end

function (op::Cheb2Coeffs2ValsOp{TR,TP})(
    coeffs::VT
) where {TR<:AbstractFloat,VT<:AbstractVector{TR},TP<:Plan}
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
