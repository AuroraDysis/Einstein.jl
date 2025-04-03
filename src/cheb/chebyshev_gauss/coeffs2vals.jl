"""
    cheb_gauss_coeffs2vals(coeffs::AbstractVector{TFC}) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    cheb_gauss_coeffs2vals_context([TF=Float64], n::Integer)(coeffs::AbstractVector{TFC})

Convert Chebyshev coefficients to values at Chebyshev points of the 1st kind.

# Performance Guide
For best performance, especially in loops or repeated calls:
```julia
ctx = cheb_gauss_coeffs2vals_context(Float64, n)
values = cheb_gauss_coeffs2vals!(ctx, coeffs)
```

# References
- [chebfun/@chebtech1/coeffs2vals.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech1/coeffs2vals.m)
"""
struct ChebyshevGaussCoeffs2ValsContext{TF<:AbstractFloat,TI<:Integer,TP<:Plan{Complex{TF}}}
    n::TI
    weights::Vector{Complex{TF}}
    tmp::Vector{Complex{TF}}
    complex_output::Vector{Complex{TF}}
    real_output::Vector{TF}
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

function _cheb_gauss_coeffs2vals_impl!(
    ctx::ChebyshevGaussCoeffs2ValsContext{TF,TI,TP}, coeffs::AbstractVector{TFC}
) where {TF<:AbstractFloat,TI<:Integer,TFC<:Union{TF,Complex{TF}},TP<:Plan{Complex{TF}}}
    (; n, weights, tmp, complex_output, fft_plan) = ctx

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
    complex_output .= @view(tmp[n:-1:1])

    # Enforce symmetry
    if isEven
        half = one(TF) / 2
        @.. complex_output = half * (complex_output + @view(complex_output[end:-1:1]))
    elseif isOdd
        half = one(TF) / 2
        @.. complex_output = half * (complex_output - @view(complex_output[end:-1:1]))
    end

    return nothing
end

function cheb_gauss_coeffs2vals!(
    ctx::ChebyshevGaussCoeffs2ValsContext{TF,TI,TP}, coeffs::AbstractVector{TF}
) where {TF<:AbstractFloat,TI<:Integer,TP<:Plan{Complex{TF}}}
    (; complex_output, real_output) = ctx
    _cheb_gauss_coeffs2vals_impl!(ctx, coeffs)
    @.. real_output = real(complex_output)
    return real_output
end

function cheb_gauss_coeffs2vals!(
    ctx::ChebyshevGaussCoeffs2ValsContext{TF,TI,TP}, coeffs::AbstractVector{Complex{TF}}
) where {TF<:AbstractFloat,TI<:Integer,TP<:Plan{Complex{TF}}}
    _cheb_gauss_coeffs2vals_impl!(ctx, coeffs)
    return ctx.complex_output
end

function cheb_gauss_coeffs2vals_context(
    ::Type{TF}, n::TI
) where {TF<:AbstractFloat,TI<:Integer}
    @argcheck n > 1 "n must be greater than 1"

    weights = Vector{Complex{TF}}(undef, 2n)
    tmp = Vector{Complex{TF}}(undef, 2n)
    complex_output = Vector{Complex{TF}}(undef, n)
    real_output = Vector{TF}(undef, n)
    fft_plan = plan_fft_measure!(tmp)

    _cheb_gauss_coeffs2vals_weights!(weights, n)

    return ChebyshevGaussCoeffs2ValsContext{TF,TI,typeof(fft_plan)}(
        n, weights, tmp, complex_output, real_output, fft_plan
    )
end

function cheb_gauss_coeffs2vals(
    coeffs::AbstractVector{TFC}
) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    n = length(coeffs)

    if n == 1
        return deepcopy(coeffs)
    end

    ctx = cheb_gauss_coeffs2vals_context(real(TFC), n)
    return cheb_gauss_coeffs2vals!(ctx, coeffs)
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
    ctx = cheb_gauss_coeffs2vals_context(TF, n)
    @inbounds for i in 1:n
        S[:, i] .= cheb_gauss_coeffs2vals!(ctx, OneElement(one(TF), i, n))
    end
    return S
end

function cheb_gauss_coeffs2vals_matrix(n::Integer)
    return cheb_gauss_coeffs2vals_matrix(Float64, n)
end

export cheb_gauss_coeffs2vals_context,
    cheb_gauss_coeffs2vals, cheb_gauss_coeffs2vals!, cheb_gauss_coeffs2vals_matrix
