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
    TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}},TI<:Integer,TP<:Plan{Complex{TF}}
}
    n::TI
    tmp::Vector{Complex{TF}}
    output::Vector{TFC}
    ifft_plan::TP
end

function _cheb_lobatto_vals2coeffs_impl!(
    ctx::ChebyshevLobattoVals2CoeffsContext{TF,TFC,TI,TP}, vals::AbstractVector{TFC}
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}},TI<:Integer,TP<:Plan{Complex{TF}}}
    (; n, tmp, output, ifft_plan) = ctx

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
    if TFC <: Real
        output[1] = real(tmp[1])
        @.. output[2:(n - 1)] = 2 * real(tmp[2:(n - 1)])
        output[n] = real(tmp[n])
    else
        output[1] = tmp[1]
        @.. output[2:(n - 1)] = 2 * tmp[2:(n - 1)]
        output[n] = tmp[n]
    end

    # Enforce exact symmetries
    if is_even
        output[2:2:end] .= 0
    elseif is_odd
        output[1:2:end] .= 0
    end

    return nothing
end

function cheb_lobatto_vals2coeffs!(
    ctx::ChebyshevLobattoVals2CoeffsContext{TF,TFC,TI,TP}, vals::AbstractVector{TFC}
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}},TI<:Integer,TP<:Plan{Complex{TF}}}
    _cheb_lobatto_vals2coeffs_impl!(ctx, vals)
    return ctx.output
end

function cheb_lobatto_vals2coeffs_context(
    ::Type{TFC}, n::TI
) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}},TI<:Integer}
    @argcheck n > 1 "n must be greater than 1"

    TF = real(TFC)
    tmp = zeros(Complex{TF}, 2n - 2)
    output = zeros(TFC, n)
    ifft_plan = plan_ifft_measure!(tmp)
    return ChebyshevLobattoVals2CoeffsContext{TF,TFC,TI,typeof(ifft_plan)}(
        n, tmp, output, ifft_plan
    )
end

function cheb_lobatto_vals2coeffs(
    vals::AbstractVector{TFC}
) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    n = length(vals)
    ctx = cheb_lobatto_vals2coeffs_context(TFC, n)
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
