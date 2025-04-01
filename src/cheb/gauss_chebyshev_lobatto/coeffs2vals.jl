"""
    coeffs2vals(coeffs::AbstractVector{TFC}) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    coeffs2vals([TF=Float64], n::Integer)(coeffs::AbstractVector{TFC})

Convert Chebyshev coefficients to values at Chebyshev points of the 2nd kind.

# Performance Guide
For best performance, especially in loops or repeated calls:
```julia
op = coeffs2vals(Float64, n)
values = op(coeffs)
```

# References
- [chebfun/@chebtech2/coeffs2vals.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/coeffs2vals.m)
"""
struct Coeffs2ValsCache{TF<:AbstractFloat}
    tmp::Vector{Complex{TF}}
    vals::Vector{Complex{TF}}
    real_vals::Vector{TF}
    fft_plan::Plan{Complex{TF}}

    function Coeffs2ValsCache{TF}(n::Integer) where {TF<:AbstractFloat}
        tmp = zeros(Complex{TF}, 2n - 2)
        vals = zeros(Complex{TF}, n)
        real_vals = zeros(TF, n)
        fft_plan = plan_fft_measure!(tmp)
        return new{TF}(tmp, vals, real_vals, fft_plan)
    end
end

function (op::Coeffs2ValsCache{TF})(
    coeffs::AbstractVector{TFC}
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}}
    type_is_float = TFC <: AbstractFloat

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

    half = one(TF) / 2
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

function coeffs2vals(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    return Coeffs2ValsCache{TF}(n)
end

function coeffs2vals(coeffs::AbstractVector{TFC}) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    n = length(coeffs)

    if n <= 1
        return deepcopy(coeffs)
    end

    op = Coeffs2ValsCache{real(TFC)}(n)
    return op(coeffs)
end

"""
    coeffs2vals_matrix([TF=Float64], n::Integer) where {TF<:AbstractFloat}

Construct the synthesis matrix S that transforms Chebyshev coefficients to function values at Chebyshev points of the 2nd kind.

# Arguments
- `TF`: Element type (defaults to Float64)
- `n`: Number of points/coefficients
"""
function coeffs2vals_matrix(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    S = Array{TF,2}(undef, n, n)
    op = Coeffs2ValsCache{TF}(n)
    @inbounds for i in 1:n
        S[:, i] = op(OneElement(one(TF), i, n))
    end
    return S
end

function coeffs2vals_matrix(n::Integer)
    return coeffs2vals_matrix(Float64, n)
end

export coeffs2vals, coeffs2vals_matrix
