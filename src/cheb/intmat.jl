"""
    cheb2_coeffs_intmat(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    cheb2_coeffs_intmat(::Type{TR}, n::TI, x_min::TR, x_max::TR) where {TR<:AbstractFloat,TI<:Integer}

Generate the Chebyshev coefficient integration matrix for spectral integration.

# Arguments
- `TR`: Type parameter for the matrix elements (e.g., Float64)
- `n`: Size of the matrix (n×n)
- `x_min`: (Optional) Lower bound of the integration interval
- `x_max`: (Optional) Upper bound of the integration interval

# Returns
- `Matrix{TR}`: The integration matrix B (n×n)

# Mathematical Background
The integration matrix B operates on Chebyshev spectral coefficients to compute
the coefficients of the indefinite integral. For a function expressed in the 
Chebyshev basis:

```math
f(x) = \\sum_{k=0}^{N-1} a_k T_k(x)
```

The indefinite integral's coefficients ``b_k`` in:

```math
\\int f(x)\\,dx = \\sum_{k=0}^{N-1} b_k T_k(x) + C
```

are computed using the matrix B: ``b = Ba``

The matrix elements are derived from the integration relation of Chebyshev polynomials:
```math
\\int T_n(x)\\,dx = \\frac{1}{2}\\left(\\frac{T_{n+1}(x)}{n+1} - \\frac{T_{n-1}(x)}{n-1}\\right)
```

When `x_min` and `x_max` are provided, the matrix is scaled for integration over [x_min, x_max].

# Examples
```julia
# Generate 8×8 integration matrix for [-1,1]
B = cheb2_coeffs_intmat(Float64, 8)

# Get Chebyshev coefficients of sin(x) using cheb2_asmat
A, _ = cheb2_asmat(Float64, 8)
x = cheb_grid(Float64, 8)
f = sin.(x)
a = A * f  # Chebyshev coefficients of sin(x)

# Compute coefficients of indefinite integral
b = B * a  # Chebyshev coefficients of -cos(x) + C
```

See also: [`cheb_grid`](@ref), [`cheb2_asmat`](@ref)
"""
function cheb2_coeffs_intmat(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    nm1 = n - 1

    B = zeros(TR, n, n)

    @inbounds begin
        B[1, 1] = 1
        B[1, 2] = -one(TR) / 4
        B[2, 1] = 1
        B[2, 3] = -one(TR) / 2
        B[3, 2] = one(TR) / 4
    end

    @inbounds for i in 3:nm1
        # upper diagonal
        B[i, i + 1] = -one(TR) / (2 * (i - 1))

        # lower diagonal
        B[i + 1, i] = one(TR) / (2 * i)

        # first row
        c = i % 2 == 0 ? -1 : 1
        B[1, i] = c * (B[i - 1, i] + B[i + 1, i])
    end

    # fix B[end-1, end]
    @inbounds B[end - 1, end] += one(TR) / (2 * n)
    @inbounds B[1, end] = (nm1 % 2 == 0 ? 1 : -1) * B[end - 1, end]

    return B
end

# Second method documentation is inherited from the main docstring
function cheb2_coeffs_intmat(
    ::Type{TR}, n::TI, x_min::TR, x_max::TR
) where {TR<:AbstractFloat,TI<:Integer}
    B = cheb2_coeffs_intmat(TR, n)
    B .*= (x_max - x_min) / 2
    return B
end

"""
    cheb2_intmat(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    cheb2_intmat(::Type{TR}, n::TI, x_min::TR, x_max::TR) where {TR<:AbstractFloat,TI<:Integer}

Generate the Chebyshev integration matrix that operates directly on function values.

# Arguments
- `TR`: Type parameter for the matrix elements (e.g., Float64)
- `n`: Size of the matrix (n×n)
- `x_min`: (Optional) Lower bound of the integration interval
- `x_max`: (Optional) Upper bound of the integration interval

# Returns
- `Matrix{TR}`: The integration matrix that operates on function values

# Mathematical Background
This matrix directly computes the indefinite integral of a function from its values
at Chebyshev points. For a function ``f(x)``, the integral is computed as:

```math
\\int f(x)\\,dx = \\mathbf{I}f
```

where ``\\mathbf{I} = S B A`` is the integration matrix composed of:
- ``A``: Analysis matrix (transform to spectral coefficients)
- ``B``: Coefficient integration matrix
- ``S``: Synthesis matrix (transform back to physical space)

This composition allows integration in physical space through:
1. Transform to spectral space (A)
2. Integrate coefficients (B)
3. Transform back to physical space (S)

# Examples
```julia
# Generate 8×8 integration matrix for [-1,1]
I = cheb2_intmat(Float64, 8)

# Get function values at Chebyshev points
x = cheb_grid(Float64, 8)
f = sin.(x)

# Compute indefinite integral (-cos(x) + C)
F = I * f

# Integration matrix for [0,π]
I_scaled = cheb2_intmat(Float64, 8, 0.0, π)
```

See also: [`cheb2_coeffs_intmat`](@ref), [`cheb2_asmat`](@ref), [`cheb_grid`](@ref)
"""
function cheb2_intmat(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    A, S = cheb2_asmat(TR, n)
    B = cheb2_coeffs_intmat(TR, n)
    return S * B * A
end

# Second method documentation is inherited from the main docstring
function cheb2_intmat(
    ::Type{TR}, n::TI, x_min::TR, x_max::TR
) where {TR<:AbstractFloat,TI<:Integer}
    A, S = cheb2_asmat(TR, n)
    B = cheb2_coeffs_intmat(TR, n, x_min, x_max)
    return S * B * A
end

export cheb2_coeffs_intmat, cheb2_intmat
