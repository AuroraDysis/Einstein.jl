"""
    gauss_chebyshev_lobatto_vals2coeffs(vals::AbstractVector{TF}) where {TF<:AbstractFloat}
    gauss_chebyshev_lobatto_vals2coeffs([TF=Float64], n::Integer)(vals::AbstractVector{TF}) where {TF<:AbstractFloat}

Convert values at Chebyshev points of the 2nd kind into Chebyshev coefficients.

# Performance Guide
For best performance, especially in loops or repeated calls:
```julia
op = gauss_chebyshev_lobatto_vals2coeffs(Float64, n)
values = op(coeffs)
```

# References
- [chebfun/@chebtech2/vals2coeffs.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/vals2coeffs.m)
"""
struct GaussChebyshevLobattoVals2CoeffsCache{TF<:AbstractFloat,TPlan<:Plan{Complex{TF}}}
    n::Integer
    tmp::Vector{Complex{TF}}
    complex_coeffs::Vector{Complex{TF}}
    real_coeffs::Vector{TF}
    ifft_plan::TPlan

    function GaussChebyshevLobattoVals2CoeffsCache{TF}(n::Integer) where {TF<:AbstractFloat}
        @argcheck n > 1 "n must be greater than 1"

        tmp = zeros(Complex{TF}, 2n - 2)
        complex_coeffs = zeros(Complex{TF}, n)
        real_coeffs = zeros(TF, n)
        ifft_plan = plan_ifft_measure!(tmp)
        return new{TF,typeof(ifft_plan)}(n, tmp, complex_coeffs, real_coeffs, ifft_plan)
    end
end

function _compute_gauss_chebyshev_lobatto_vals2coeffs!(
    op::GaussChebyshevLobattoVals2CoeffsCache{TF}, vals::AbstractVector{TFC}
) where {TF<:AbstractFloat,TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    (; n, tmp, complex_coeffs, ifft_plan) = op

    @argcheck length(vals) == n "vals must have length n"

    # Determine if vals are even or odd symmetric with tolerance
    tol = sqrt(eps(TF))
    is_even = all(i -> abs(vals[i] - vals[n - i + 1]) < tol, 1:(n ÷ 2))
    is_odd = all(i -> abs(vals[i] + vals[n - i + 1]) < tol, 1:(n ÷ 2))

    # Mirror the values
    @inbounds for i in 1:(n - 1)
        tmp[i] = vals[n - i + 1]  # descending part
        tmp[n - 1 + i] = vals[i]  # ascending part
    end

    # Perform inverse FFT on the mirrored data
    ifft_plan * tmp

    # Scale the interior coefficients
    @inbounds begin
        complex_coeffs[1] = tmp[1]
        @. complex_coeffs[2:(n - 1)] = 2 * @view(tmp[2:(n - 1)])
        complex_coeffs[n] = tmp[n]
    end

    # Enforce exact symmetries
    if is_even
        complex_coeffs[2:2:end] .= 0
    elseif is_odd
        complex_coeffs[1:2:end] .= 0
    end

    return nothing
end

function (op::GaussChebyshevLobattoVals2CoeffsCache{TF})(
    vals::AbstractVector{TF}
) where {TF<:AbstractFloat}
    (; complex_coeffs, real_coeffs) = op
    _compute_gauss_chebyshev_lobatto_vals2coeffs!(op, vals)
    @. real_coeffs = real(complex_coeffs)
    return real_coeffs
end

function (op::GaussChebyshevLobattoVals2CoeffsCache{TF})(
    vals::AbstractVector{Complex{TF}}
) where {TF<:AbstractFloat}
    (; complex_coeffs) = op
    _compute_gauss_chebyshev_lobatto_vals2coeffs!(op, vals)
    return complex_coeffs
end

function gauss_chebyshev_lobatto_vals2coeffs(
    ::Type{TF}, n::Integer
) where {TF<:AbstractFloat}
    return GaussChebyshevLobattoVals2CoeffsCache{TF}(n)
end

function gauss_chebyshev_lobatto_vals2coeffs(
    vals::AbstractVector{TFC}
) where {TFC<:Union{AbstractFloat,Complex{<:AbstractFloat}}}
    n = length(vals)
    op = GaussChebyshevLobattoVals2CoeffsCache{real(TFC)}(n)
    return op(vals)
end

"""
    gauss_chebyshev_lobatto_vals2coeffs_matrix([TF=Float64], n::Integer) where {TF<:AbstractFloat}

Construct the analysis matrix A that transforms function values at Chebyshev points of the 2nd kind to Chebyshev coefficients.

# Arguments
- `TF`: Element type (defaults to Float64)
- `n`: Number of points/coefficients
"""
function gauss_chebyshev_lobatto_vals2coeffs_matrix(
    ::Type{TF}, n::Integer
) where {TF<:AbstractFloat}
    @argcheck n > 1 "n must be greater than 1"

    A = Array{TF,2}(undef, n, n)
    op = GaussChebyshevLobattoVals2CoeffsCache{TF}(n)
    @inbounds for i in 1:n
        A[:, i] = op(OneElement(one(TF), i, n))
    end
    return A
end

function gauss_chebyshev_lobatto_vals2coeffs_matrix(n::Integer)
    return gauss_chebyshev_lobatto_vals2coeffs_matrix(Float64, n)
end

export gauss_chebyshev_lobatto_vals2coeffs, gauss_chebyshev_lobatto_vals2coeffs_matrix
