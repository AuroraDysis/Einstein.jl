"""
    cheb_lobatto_vals2coeffs(vals::AbstractVector{TF}) where {TF<:AbstractFloat}
    cheb_lobatto_vals2coeffs_context([TF=Float64], n::Integer)(vals::AbstractVector{TF}) where {TF<:AbstractFloat}

Convert values at Chebyshev points of the 2nd kind into Chebyshev coefficients.

# Performance Guide
For best performance, especially in loops or repeated calls:
```julia
ctx = cheb_lobatto_vals2coeffs_context(Float64, n)
coeffs = cheb_lobatto_vals2coeffs!(ctx, values)
```

# References
- [chebfun/@chebtech2/vals2coeffs.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/vals2coeffs.m)
"""
struct ChebyshevLobattoVals2CoeffsContext{
    TF<:AbstractFloat,TI<:Integer,TP<:Plan{Complex{TF}}
}
    n::TI
    tmp::Vector{Complex{TF}}
    complex_output::Vector{Complex{TF}}
    real_output::Vector{TF}
    ifft_plan::TP
end

function _cheb_lobatto_vals2coeffs_impl!(
    ctx::ChebyshevLobattoVals2CoeffsContext{TF,TI,TP}, vals::AbstractVector{TFC}
) where {TF<:AbstractFloat,TI<:Integer,TFC<:Union{TF,Complex{TF}},TP<:Plan{Complex{TF}}}
    (; n, tmp, complex_output, ifft_plan) = ctx

    @argcheck length(vals) == n "vals must have length n"

    # Determine if vals are even or odd symmetric with tolerance
    tol = sqrt(eps(TF))
    is_even = all(i -> abs(vals[i] - vals[n - i + 1]) < tol, 1:(n ÷ 2))
    is_odd = all(i -> abs(vals[i] + vals[n - i + 1]) < tol, 1:(n ÷ 2))

    # Mirror the values
    tmp[1:(n - 1)] .= @view(vals[n:-1:2])
    tmp[n:(2n - 2)] .= @view(vals[1:(n - 1)])

    # Perform inverse FFT on the mirrored data
    ifft_plan * tmp

    # Scale the interior coefficients
    complex_output[1] = tmp[1]
    @.. complex_output[2:(n - 1)] = 2 * @view(tmp[2:(n - 1)])
    complex_output[n] = tmp[n]

    # Enforce exact symmetries
    if is_even
        complex_output[2:2:end] .= 0
    elseif is_odd
        complex_output[1:2:end] .= 0
    end

    return nothing
end

function cheb_lobatto_vals2coeffs!(
    ctx::ChebyshevLobattoVals2CoeffsContext{TF,TI,TP}, vals::AbstractVector{TF}
) where {TF<:AbstractFloat,TI<:Integer,TP<:Plan{Complex{TF}}}
    (; complex_output, real_output) = ctx
    _cheb_lobatto_vals2coeffs_impl!(ctx, vals)
    @. real_output = real(complex_output)
    return real_output
end

function cheb_lobatto_vals2coeffs!(
    ctx::ChebyshevLobattoVals2CoeffsContext{TF,TI,TP}, vals::AbstractVector{Complex{TF}}
) where {TF<:AbstractFloat,TI<:Integer,TP<:Plan{Complex{TF}}}
    (; complex_output) = ctx
    _cheb_lobatto_vals2coeffs_impl!(ctx, vals)
    return complex_output
end

function cheb_lobatto_vals2coeffs_context(
    ::Type{TF}, n::TI
) where {TF<:AbstractFloat,TI<:Integer}
    @argcheck n > 1 "n must be greater than 1"

    tmp = zeros(Complex{TF}, 2n - 2)
    complex_output = zeros(Complex{TF}, n)
    real_output = zeros(TF, n)
    ifft_plan = plan_ifft_measure!(tmp)
    return ChebyshevLobattoVals2CoeffsContext{TF,TI,typeof(ifft_plan)}(
        n, tmp, complex_output, real_output, ifft_plan
    )
end

function cheb_lobatto_vals2coeffs(
    vals::AbstractVector{TFC}
) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    n = length(vals)
    ctx = cheb_lobatto_vals2coeffs_context(real(TFC), n)
    return cheb_lobatto_vals2coeffs!(ctx, vals)
end

"""
    cheb_lobatto_vals2coeffs_matrix([TF=Float64], n::Integer) where {TF<:AbstractFloat}

Construct the analysis matrix A that transforms function values at Chebyshev points of the 2nd kind to Chebyshev coefficients.

# Arguments
- `TF`: Element type (defaults to Float64)
- `n`: Number of points/coefficients
"""
function cheb_lobatto_vals2coeffs_matrix(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    @argcheck n > 0 "n must be greater than 0"

    if n == 1
        return ones(TF, 1, 1)
    end

    A = Array{TF,2}(undef, n, n)
    ctx = cheb_lobatto_vals2coeffs_context(TF, n)
    @inbounds for i in 1:n
        A[:, i] .= cheb_lobatto_vals2coeffs!(ctx, OneElement(one(TF), i, n))
    end
    return A
end

function cheb_lobatto_vals2coeffs_matrix(n::Integer)
    return cheb_lobatto_vals2coeffs_matrix(Float64, n)
end

export cheb_lobatto_vals2coeffs_context,
    cheb_lobatto_vals2coeffs, cheb_lobatto_vals2coeffs!, cheb_lobatto_vals2coeffs_matrix
