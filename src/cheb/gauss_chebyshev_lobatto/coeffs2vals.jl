"""
    gauss_chebyshev_lobatto_coeffs2vals(coeffs::AbstractVector{TFC}) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    gauss_chebyshev_lobatto_coeffs2vals([TF=Float64], n::Integer)(coeffs::AbstractVector{TFC})

Convert Chebyshev coefficients to values at Chebyshev points of the 2nd kind.

# Performance Guide
For best performance, especially in loops or repeated calls:
```julia
op = gauss_chebyshev_lobatto_coeffs2vals(Float64, n)
values = op(coeffs)
```

# References
- [chebfun/@chebtech2/coeffs2vals.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/coeffs2vals.m)
"""
struct GaussChebyshevLobattoCoeffs2ValsCache{TF<:AbstractFloat,TPlan<:Plan{Complex{TF}}}
    n::Integer
    tmp::Vector{Complex{TF}}
    complex_output::Vector{Complex{TF}}
    real_output::Vector{TF}
    fft_plan::TPlan

    function GaussChebyshevLobattoCoeffs2ValsCache{TF}(n::Integer) where {TF<:AbstractFloat}
        @argcheck n > 1 "n must be greater than 1"

        tmp = zeros(Complex{TF}, 2n - 2)
        complex_output = zeros(Complex{TF}, n)
        real_output = zeros(TF, n)
        fft_plan = plan_fft_measure!(tmp)
        return new{TF,typeof(fft_plan)}(n, tmp, complex_output, real_output, fft_plan)
    end
end

function _compute_gauss_chebyshev_lobatto_coeffs2vals!(
    op::GaussChebyshevLobattoCoeffs2ValsCache{TF}, coeffs::AbstractVector{TFC}
) where {TF<:AbstractFloat,TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    (; n, tmp, complex_output, fft_plan) = op

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
        complex_output .= @view(tmp[n:-1:1])
    end

    # In-place symmetry enforcement
    if isEven
        half = one(TF) / 2
        @. complex_output = half * (complex_output + @view(complex_output[n:-1:1]))
    elseif isOdd
        half = one(TF) / 2
        @. complex_output = half * (complex_output - @view(complex_output[n:-1:1]))
    end

    return nothing
end

function (op::GaussChebyshevLobattoCoeffs2ValsCache{TF})(
    coeffs::AbstractVector{TF}
) where {TF<:AbstractFloat}
    (; complex_output, real_output) = op
    _compute_gauss_chebyshev_lobatto_coeffs2vals!(op, coeffs)
    @. real_output = real(complex_output)
    return real_output
end

function (op::GaussChebyshevLobattoCoeffs2ValsCache{TF})(
    coeffs::AbstractVector{Complex{TF}}
) where {TF<:AbstractFloat}
    _compute_gauss_chebyshev_lobatto_coeffs2vals!(op, coeffs)
    return op.complex_output
end

function gauss_chebyshev_lobatto_coeffs2vals(
    ::Type{TF}, n::Integer
) where {TF<:AbstractFloat}
    return GaussChebyshevLobattoCoeffs2ValsCache{TF}(n)
end

function gauss_chebyshev_lobatto_coeffs2vals(
    coeffs::AbstractVector{TFC}
) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    n = length(coeffs)
    op = GaussChebyshevLobattoCoeffs2ValsCache{real(TFC)}(n)
    return op(coeffs)
end

"""
    gauss_chebyshev_lobatto_coeffs2vals_matrix([TF=Float64], n::Integer) where {TF<:AbstractFloat}

Construct the synthesis matrix S that transforms Chebyshev coefficients to function values at Chebyshev points of the 2nd kind.

# Arguments
- `TF`: Element type (defaults to Float64)
- `n`: Number of points/coefficients
"""
function gauss_chebyshev_lobatto_coeffs2vals_matrix(
    ::Type{TF}, n::Integer
) where {TF<:AbstractFloat}
    S = Array{TF,2}(undef, n, n)
    op = GaussChebyshevLobattoCoeffs2ValsCache{TF}(n)
    @inbounds for i in 1:n
        S[:, i] = op(OneElement(one(TF), i, n))
    end
    return S
end

function gauss_chebyshev_lobatto_coeffs2vals_matrix(n::Integer)
    return gauss_chebyshev_lobatto_coeffs2vals_matrix(Float64, n)
end

export gauss_chebyshev_lobatto_coeffs2vals, gauss_chebyshev_lobatto_coeffs2vals_matrix
