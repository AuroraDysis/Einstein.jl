"""
    cheb_gauss_vals2coeffs(vals::AbstractVector{TF}) where {TF<:AbstractFloat}
    cheb_gauss_vals2coeffs_context([TF=Float64], n::Integer)(vals::AbstractVector{TF}) where {TF<:AbstractFloat}

Convert values at Chebyshev points of the 1st kind into Chebyshev coefficients.

# Performance Guide
For best performance, especially in loops or repeated calls:
```julia
ctx = cheb_gauss_vals2coeffs_context(Float64, n)
coeffs = cheb_gauss_vals2coeffs!(ctx, values)
```

# References
- [chebfun/@chebtech1/vals2coeffs.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech1/vals2coeffs.m)
"""
struct ChebyshevGaussVals2CoeffsContext{TF<:AbstractFloat,TI<:Integer,TP<:Plan{Complex{TF}}}
    n::TI
    weights::Vector{Complex{TF}}
    tmp::Vector{Complex{TF}}
    complex_output::Vector{Complex{TF}}
    real_output::Vector{TF}
    ifft_plan::TP
end

function _cheb_gauss_vals2coeffs_weights!(
    weights::AbstractVector{Complex{TF}}, n::Integer
) where {TF<:AbstractFloat}
    @inbounds begin
        im_pi_over_2n = im * convert(TF, π) / (2n)
        @.. weights = 2 * exp((0:(n - 1)) * im_pi_over_2n)
        weights[1] /= 2
    end

    return nothing
end

function _cheb_gauss_vals2coeffs_impl!(
    ctx::ChebyshevGaussVals2CoeffsContext{TF,TI,TP}, vals::AbstractVector{TFC}
) where {TF<:AbstractFloat,TI<:Integer,TFC<:Union{TF,Complex{TF}},TP<:Plan{Complex{TF}}}
    (; n, weights, tmp, complex_output, ifft_plan) = ctx

    @argcheck length(vals) == n "vals must have length n"

    # Check for symmetry with tolerance
    tol = sqrt(eps(TF))
    isEven = all(i -> abs(vals[i] - vals[n - i + 1]) < tol, 1:(n ÷ 2))
    isOdd = all(i -> abs(vals[i] + vals[n - i + 1]) < tol, 1:(n ÷ 2))

    # Mirror the values for FFT
    tmp[1:n] .= @view(vals[n:-1:1])
    tmp[(n + 1):(2n)] .= vals

    ifft_plan * tmp

    @.. complex_output = @view(tmp[1:n]) * weights

    # adjust coefficients for symmetry
    if isEven
        complex_output[2:2:end] .= 0
    elseif isOdd
        complex_output[1:2:end] .= 0
    end

    return nothing
end

function cheb_gauss_vals2coeffs!(
    ctx::ChebyshevGaussVals2CoeffsContext{TF,TI,TP}, values::AbstractVector{TF}
) where {TF<:AbstractFloat,TI<:Integer,TP<:Plan{Complex{TF}}}
    (; complex_output, real_output) = ctx
    _cheb_gauss_vals2coeffs_impl!(ctx, values)
    @. real_output = real(complex_output)
    return real_output
end

function cheb_gauss_vals2coeffs!(
    ctx::ChebyshevGaussVals2CoeffsContext{TF,TI,TP}, values::AbstractVector{Complex{TF}}
) where {TF<:AbstractFloat,TI<:Integer,TP<:Plan{Complex{TF}}}
    (; complex_output) = ctx
    _cheb_gauss_vals2coeffs_impl!(ctx, values)
    return complex_output
end

function cheb_gauss_vals2coeffs_context(
    ::Type{TF}, n::TI
) where {TF<:AbstractFloat,TI<:Integer}
    @argcheck n > 1 "n must be greater than 1"

    weights = Vector{Complex{TF}}(undef, n)
    tmp = Vector{Complex{TF}}(undef, 2n)
    complex_output = Vector{Complex{TF}}(undef, n)
    real_output = Vector{TF}(undef, n)
    ifft_plan = plan_ifft_measure!(tmp)

    _cheb_gauss_vals2coeffs_weights!(weights, n)

    return ChebyshevGaussVals2CoeffsContext{TF,TI,typeof(ifft_plan)}(
        n, weights, tmp, complex_output, real_output, ifft_plan
    )
end

function cheb_gauss_vals2coeffs(
    values::AbstractVector{TFC}
) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    n = length(values)
    ctx = cheb_gauss_vals2coeffs_context(real(TFC), n)
    return cheb_gauss_vals2coeffs!(ctx, values)
end

"""
    cheb_gauss_vals2coeffs_matrix([TF=Float64], n::Integer) where {TF<:AbstractFloat}

Construct the analysis matrix A that transforms function values at Chebyshev points of the 1st kind to Chebyshev coefficients.

# Arguments
- `TF`: Element type (defaults to Float64)
- `n`: Number of points/coefficients
"""
function cheb_gauss_vals2coeffs_matrix(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    @argcheck n > 0 "n must be greater than 0"

    if n == 1
        return ones(TF, 1, 1)
    end

    A = Array{TF,2}(undef, n, n)
    ctx = cheb_gauss_vals2coeffs_context(TF, n)
    @inbounds for i in 1:n
        A[:, i] .= cheb_gauss_vals2coeffs!(ctx, OneElement(one(TF), i, n))
    end
    return A
end

function cheb_gauss_vals2coeffs_matrix(n::Integer)
    return cheb_gauss_vals2coeffs_matrix(Float64, n)
end

export cheb_gauss_vals2coeffs_context,
    cheb_gauss_vals2coeffs, cheb_gauss_vals2coeffs!, cheb_gauss_vals2coeffs_matrix
