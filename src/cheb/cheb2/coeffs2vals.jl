"""
    cheb2_coeffs2vals(coeffs::AbstractVector{TR}) where {TR<:AbstractFloat}
    op::Cheb2Coeffs2ValsOp{[TR=Float64]}(n::Integer)(coeffs::AbstractVector{TR}) where {TR<:AbstractFloat}

Convert Chebyshev coefficients to values at Chebyshev points of the 2nd kind.

# Performance Guide
For best performance, especially in loops or repeated calls:
```julia
op = Cheb2Coeffs2ValsOp{Float64}(n)
values = op(coeffs)
```

# References
- [chebfun/@chebtech2/coeffs2vals.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/coeffs2vals.m)
"""
struct Cheb2Coeffs2ValsOp{TR<:AbstractFloat}
    tmp::Vector{Complex{TR}}
    vals::Vector{Complex{TR}}
    real_vals::Vector{TR}
    fft_plan::Plan{Complex{TR}}

    function Cheb2Coeffs2ValsOp{TR}(n::Integer) where {TR<:AbstractFloat}
        tmp = zeros(Complex{TR}, 2n - 2)
        vals = zeros(Complex{TR}, n)
        real_vals = zeros(TR, n)
        fft_plan = plan_fft_measure!(tmp)
        return new{TR}(tmp, vals, real_vals, fft_plan)
    end

    function Cheb2Coeffs2ValsOp(n::Integer)
        return Cheb2Coeffs2ValsOp{Float64}(n)
    end
end

function (op::Cheb2Coeffs2ValsOp{TR})(
    coeffs::AbstractVector{TR}
) where {TR<:AbstractFloatOrComplex}
    type_is_float = typeisfloat(TR)

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
            vals[i] = tmp[n - i + 1]
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

    if type_is_float
        @inbounds for k in 1:n
            op.real_vals[k] = real(op.vals[k])
        end
        return op.real_vals
    else
        return op.vals
    end
end

function cheb2_coeffs2vals(coeffs::AbstractVector{TR}) where {TR<:AbstractFloatOrComplex}
    n = length(coeffs)

    if n <= 1
        return deepcopy(coeffs)
    end

    op = Cheb2Coeffs2ValsOp{real(TR)}(n)
    return op(coeffs)
end

export cheb2_coeffs2vals, Cheb2Coeffs2ValsOp
