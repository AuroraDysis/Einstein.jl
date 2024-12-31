"""
    cheb2_asmat(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}

Generate the analysis and synthesis matrices for Chebyshev spectral methods.

# Arguments
- `TR`: Type parameter for the matrix elements (e.g., Float64)
- `n`: Size of the matrices (n×n)

# Returns
- `Tuple{Matrix{TR}, Matrix{TR}}`: A tuple containing:
  - Analysis matrix A (n×n)
  - Synthesis matrix S (n×n)

# Mathematical Background
The analysis and synthesis matrices are used for transforming between physical and
spectral spaces in Chebyshev spectral methods.

For a function ``f(x)`` evaluated at Chebyshev points, these matrices allow:
- Transformation to spectral coefficients: ``\\hat{f} = Af``
- Transformation back to physical space: ``f = S\\hat{f}``

The matrices are constructed using:
```math
S_{ij} = \\epsilon_j \\cos\\left(\\frac{\\pi i j}{N-1}\\right)
```
```math
A_{ji} = \\frac{2c_ic_j}{N-1}S_{ij}
```
where:
- ``c_k = \\begin{cases} 1/2 & k=0 \\text{ or } k=N-1 \\\\ 1 & \\text{otherwise} \\end{cases}``
- ``\\epsilon_j = (-1)^j``
- ``i,j = 0,\\ldots,N-1``

# Examples
```julia
# Generate 8×8 analysis and synthesis matrices with Float64 precision
A, S = cheb2_asmat(Float64, 8)

# Transform function values to spectral coefficients
f_values = [sin(x) for x in cheb2_grid(Float64, 8)]
f_coeffs = A * f_values

# Transform back to physical space
f_recovered = S * f_coeffs
```

See also: [`cheb2_grid`](@ref)
"""
function cheb2_asmat(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    nm1 = n - 1

    # analysis matrix
    A = zeros(TR, n, n)

    # synthesis matrix
    S = zeros(TR, n, n)

    inv_nm1 = one(TR) / nm1
    pi_nm1 = convert(TR, π) * inv_nm1

    for j in 0:nm1
        pm = if j % 2 == 0
            1
        else
            -1
        end

        for i in 0:nm1
            S[i + 1, j + 1] = pm * cos(i * j * pi_nm1)

            ci = if i == 0 || i == nm1
                one(TR) / 2
            else
                one(TR)
            end

            cj = if j == 0 || j == nm1
                one(TR) / 2
            else
                one(TR)
            end

            A[j + 1, i + 1] = 2 * ci * cj * inv_nm1 * S[i + 1, j + 1]
        end
    end

    return A, S
end

export cheb2_asmat
