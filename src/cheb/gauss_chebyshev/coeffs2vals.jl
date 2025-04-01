"""
    gauss_chebyshev_coeffs2vals(coeffs::AbstractVector{TFC}) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    gauss_chebyshev_coeffs2vals([TF=Float64], n::Integer)(coeffs::AbstractVector{TFC})

Convert Chebyshev coefficients to values at Chebyshev points of the 1st kind.

# Performance Guide
For best performance, especially in loops or repeated calls:
```julia
op = gauss_chebyshev_coeffs2vals(Float64, n)
values = op(coeffs)
```

# References
- [chebfun/@chebtech1/coeffs2vals.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech1/coeffs2vals.m)
"""
struct GaussChebyshevCoeffs2ValsCache{TF<:AbstractFloat,TPlan<:Plan{Complex{TF}}}
    weights::Vector{Complex{TF}}
    tmp::Vector{Complex{TF}}
    complex_values::Vector{Complex{TF}}
    real_values::Vector{TF}
    fft_plan::TPlan

    function GaussChebyshevCoeffs2ValsCache{TF}(n::Integer) where {TF<:AbstractFloat}
        # Precompute weights
        weights = Vector{Complex{TF}}(undef, 2n)
        @inbounds begin
            half = one(TF) / 2
            m_im_pi_over_2n = -im * convert(TF, π) / (2n)
            for k in 0:(n - 1)
                weights[k + 1] = exp(k * m_im_pi_over_2n) * half
            end
            weights[1] *= 2
            weights[n + 1] = 0
            for k in (n + 1):(2n - 1)
                weights[k + 1] = -exp(k * m_im_pi_over_2n) * half
            end
        end
        tmp = Vector{Complex{TF}}(undef, 2n)
        complex_values = Vector{Complex{TF}}(undef, n)
        real_values = Vector{TF}(undef, n)
        fft_plan = plan_fft_measure!(tmp)
        return new{TF,typeof(fft_plan)}(weights, tmp, complex_values, real_values, fft_plan)
    end
end

function (op::GaussChebyshevCoeffs2ValsCache{TF})(
    coeffs::AbstractVector{TFC}
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}}
    type_is_float = TFC <: AbstractFloat

    n = length(coeffs)
    if n <= 1
        if type_is_float
            op.real_values .= coeffs
            return op.real_values
        else
            op.complex_values .= coeffs
            return op.complex_values
        end
    end

    weights = op.weights
    tmp = op.tmp
    complex_values = op.complex_values
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

    @inbounds begin
        tmp .*= weights
        fft_plan * tmp
    end

    # Extract values
    @inbounds for i in 1:n
        complex_values[i] = tmp[n - i + 1]
    end

    # Enforce symmetry if needed
    if isEven || isOdd
        half = one(TF) / 2
        @inbounds for i in 1:div(n, 2)
            j = n - i + 1
            if isEven
                s = complex_values[i] + complex_values[j]
                complex_values[i] = half * s
                complex_values[j] = half * s
            else
                d = complex_values[i] - complex_values[j]
                complex_values[i] = half * d
                complex_values[j] = -half * d
            end
        end
    end

    if type_is_float
        @inbounds for k in 1:n
            op.real_values[k] = real(op.complex_values[k])
        end
        return op.real_values
    else
        return op.complex_values
    end
end

function gauss_chebyshev_coeffs2vals(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    return GaussChebyshevCoeffs2ValsCache{TF}(n)
end

function gauss_chebyshev_coeffs2vals(
    coeffs::AbstractVector{TFC}
) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    n = length(coeffs)

    if n <= 1
        return deepcopy(coeffs)
    end

    op = GaussChebyshevCoeffs2ValsCache{real(TFC)}(n)
    return op(coeffs)
end

"""
    gauss_chebyshev_coeffs2vals_matrix([TF=Float64], n::Integer) where {TF<:AbstractFloat}

Construct the synthesis matrix S that transforms Chebyshev coefficients to function values at Chebyshev points of the 1st kind.

# Arguments
- `TF`: Element type (defaults to Float64)
- `n`: Number of points/coefficients
"""
function gauss_chebyshev_coeffs2vals_matrix(
    ::Type{TF}, n::Integer
) where {TF<:AbstractFloat}
    S = Array{TF,2}(undef, n, n)
    op = GaussChebyshevCoeffs2ValsCache{TF}(n)
    @inbounds for i in 1:n
        S[:, i] = op(OneElement(one(TF), i, n))
    end
    return S
end

function gauss_chebyshev_coeffs2vals_matrix(n::Integer)
    return gauss_chebyshev_coeffs2vals_matrix(Float64, n)
end

export gauss_chebyshev_coeffs2vals, gauss_chebyshev_coeffs2vals_matrix
