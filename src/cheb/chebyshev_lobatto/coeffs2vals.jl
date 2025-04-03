"""
    cheb_lobatto_coeffs2vals(coeffs::AbstractVector{TFC}) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    cheb_lobatto_coeffs2vals_context([TF=Float64], n::Integer)(coeffs::AbstractVector{TFC})

Convert Chebyshev coefficients to values at Chebyshev points of the 2nd kind.

# Performance Guide
For best performance, especially in loops or repeated calls:
```julia
ctx = cheb_lobatto_coeffs2vals_context(Float64, n)
values = cheb_lobatto_coeffs2vals!(ctx, coeffs)
```

# References
- [chebfun/@chebtech2/coeffs2vals.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/coeffs2vals.m)
"""
struct ChebyshevLobattoCoeffs2ValsContext{
    TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}},TI<:Integer,TP<:Plan{Complex{TF}}
}
    n::TI
    tmp::Vector{Complex{TF}}
    output::Vector{TFC}
    fft_plan::TP
end

function _cheb_lobatto_coeffs2vals_impl!(
    ctx::ChebyshevLobattoCoeffs2ValsContext{TF,TFC,TI,TP}, coeffs::AbstractVector{TFC}
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}},TI<:Integer,TP<:Plan{Complex{TF}}}
    (; n, tmp, output, fft_plan) = ctx

    @argcheck length(coeffs) == n "coeffs must have length n"

    # Determine which columns are purely even or purely odd based on middle coefficients
    tol = sqrt(eps(TF))
    isEven = all(x -> abs(x) < tol, @view(coeffs[2:2:end]))
    isOdd = all(x -> abs(x) < tol, @view(coeffs[1:2:end]))

    half = one(TF) / 2
    @inbounds begin
        tmp[1] = coeffs[1]
        for i in 2:(n - 1)
            hc = half * coeffs[i]
            tmp[i] = hc
            tmp[2n - i] = hc
        end
        tmp[n] = coeffs[n]

        # FFT into complex_output
        fft_plan * @view(tmp[1:(2n - 2)])

        # Flip and truncate:
        if TFC <: Real
            @.. output = real(tmp[n:-1:1])
        else
            @.. output = tmp[n:-1:1]
        end
    end

    # In-place symmetry enforcement
    if isEven
        half = one(TF) / 2
        @.. output = half * (output + @view(output[n:-1:1]))
    elseif isOdd
        half = one(TF) / 2
        @.. output = half * (output - @view(output[n:-1:1]))
    end

    return nothing
end

function cheb_lobatto_coeffs2vals!(
    ctx::ChebyshevLobattoCoeffs2ValsContext{TF,TFC,TI,TP}, coeffs::AbstractVector{TFC}
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}},TI<:Integer,TP<:Plan{Complex{TF}}}
    _cheb_lobatto_coeffs2vals_impl!(ctx, coeffs)
    return ctx.output
end

function cheb_lobatto_coeffs2vals_context(
    ::Type{TFC}, n::TI
) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}},TI<:Integer}
    @argcheck n > 1 "n must be greater than 1"

    TF = real(TFC)
    tmp = zeros(Complex{TF}, 2n - 2)
    output = zeros(TFC, n)
    fft_plan = plan_fft_measure!(tmp)
    return ChebyshevLobattoCoeffs2ValsContext{TF,TFC,TI,typeof(fft_plan)}(
        n, tmp, output, fft_plan
    )
end

function cheb_lobatto_coeffs2vals(
    coeffs::AbstractVector{TFC}
) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    n = length(coeffs)
    ctx = cheb_lobatto_coeffs2vals_context(TFC, n)
    return cheb_lobatto_coeffs2vals!(ctx, coeffs)
end

"""
    cheb_lobatto_coeffs2vals_matrix([TF=Float64], n::Integer) where {TF<:AbstractFloat}

Construct the synthesis matrix S that transforms Chebyshev coefficients to function values at Chebyshev points of the 2nd kind.

# Arguments
- `TF`: Element type (defaults to Float64)
- `n`: Number of points/coefficients
"""
function cheb_lobatto_coeffs2vals_matrix(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    @argcheck n > 0 "n must be greater than 0"

    if n == 1
        return ones(TF, 1, 1)
    end

    S = Array{TF,2}(undef, n, n)
    ctx = cheb_lobatto_coeffs2vals_context(TF, n)
    @inbounds for i in 1:n
        S[:, i] .= cheb_lobatto_coeffs2vals!(ctx, OneElement(one(TF), i, n))
    end
    return S
end

function cheb_lobatto_coeffs2vals_matrix(n::Integer)
    return cheb_lobatto_coeffs2vals_matrix(Float64, n)
end

export cheb_lobatto_coeffs2vals_context,
    cheb_lobatto_coeffs2vals, cheb_lobatto_coeffs2vals!, cheb_lobatto_coeffs2vals_matrix
