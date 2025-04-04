"""
    cheb_gauss_coeffs2vals(coeffs::AbstractVector{TFC}) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    cheb_gauss_coeffs2vals_plan([TFC=Float64], n::Integer)(coeffs::AbstractVector{TFC}) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}

Convert Chebyshev coefficients to values at Chebyshev points of the 1st kind.

# Performance Guide
For best performance, especially in loops or repeated calls:
```julia
plan = cheb_gauss_coeffs2vals_plan(Float64, n)
values = plan * coeffs
```

# References
- [chebfun/@chebtech1/coeffs2vals.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech1/coeffs2vals.m)
"""
struct ChebyshevGaussCoeffs2ValsPlan{
    TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}},TI<:Integer,TP<:Plan{Complex{TF}}
}
    n::TI
    weights::Vector{Complex{TF}}
    tmp::Vector{Complex{TF}}
    output::Vector{TFC}
    fft_plan::TP
end

function _cheb_gauss_coeffs2vals_weights!(
    weights::AbstractVector{Complex{TF}}, n::Integer
) where {TF<:AbstractFloat}
    half = one(TF) / 2
    m_im_pi_over_2n = -im * convert(TF, π) / (2n)
    @.. weights[1:(2n)] = exp((0:(2n - 1)) * m_im_pi_over_2n) * half
    weights[1] *= 2
    weights[n + 1] = 0
    return weights[(n + 2):(2n)] .*= -1
end

function *(
    plan::ChebyshevGaussCoeffs2ValsPlan{TF,TFC,TI,TP}, x::AbstractVector{TFC}
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}},TI<:Integer,TP<:Plan{Complex{TF}}}
    (; n, weights, tmp, output, fft_plan) = plan

    @argcheck length(coeffs) == n "coeffs must have length n"

    # Check for symmetry
    tol = sqrt(eps(TF))
    isEven = all(x -> abs(x) < tol, @view(coeffs[2:2:end]))
    isOdd = all(x -> abs(x) < tol, @view(coeffs[1:2:end]))

    # Copy coefficients and mirror
    tmp[1:n] .= coeffs
    tmp[n + 1] = 1
    tmp[(n + 2):(2n)] .= @view(coeffs[n:-1:2])

    tmp .*= weights
    fft_plan * tmp

    # Truncate and flip the order
    if TFC <: Real
        @.. output = real(tmp[n:-1:1])
    else
        @.. output = tmp[n:-1:1]
    end

    # Enforce symmetry
    if isEven
        half = one(TF) / 2
        @.. output = half * (output + output[end:-1:1])
    elseif isOdd
        half = one(TF) / 2
        @.. output = half * (output - output[end:-1:1])
    end

    return output
end

function cheb_gauss_coeffs2vals_plan(
    ::Type{TFC}, n::TI
) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}},TI<:Integer}
    @argcheck n > 1 "n must be greater than 1"

    TF = real(TFC)
    weights = Vector{Complex{TF}}(undef, 2n)
    tmp = Vector{Complex{TF}}(undef, 2n)
    output = Vector{TFC}(undef, n)
    fft_plan = plan_fft_measure!(tmp)

    _cheb_gauss_coeffs2vals_weights!(weights, n)

    return ChebyshevGaussCoeffs2ValsPlan{TF,TFC,TI,typeof(fft_plan)}(
        n, weights, tmp, output, fft_plan
    )
end

function cheb_gauss_coeffs2vals_plan(n::TI) where {TI<:Integer}
    return cheb_gauss_coeffs2vals_plan(Float64, n)
end

function cheb_gauss_coeffs2vals(
    coeffs::AbstractVector{TFC}
) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    n = length(coeffs)

    if n == 1
        return deepcopy(coeffs)
    end

    plan = cheb_gauss_coeffs2vals_plan(TFC, n)
    return plan * coeffs
end

"""
    cheb_gauss_coeffs2vals_matrix([TF=Float64], n::Integer) where {TF<:AbstractFloat}

Construct the synthesis matrix S that transforms Chebyshev coefficients to function values at Chebyshev points of the 1st kind.

# Arguments
- `TF`: Element type (defaults to Float64)
- `n`: Number of points/coefficients
"""
function cheb_gauss_coeffs2vals_matrix(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    @argcheck n > 0 "n must be greater than 0"

    if n == 1
        return ones(TF, 1, 1)
    end

    S = Array{TF,2}(undef, n, n)
    plan = cheb_gauss_coeffs2vals_plan(TF, n)
    @inbounds for i in 1:n
        S[:, i] .= plan * OneElement(one(TF), i, n)
    end
    return S
end

function cheb_gauss_coeffs2vals_matrix(n::Integer)
    return cheb_gauss_coeffs2vals_matrix(Float64, n)
end

export cheb_gauss_coeffs2vals_plan, cheb_gauss_coeffs2vals, cheb_gauss_coeffs2vals_matrix
