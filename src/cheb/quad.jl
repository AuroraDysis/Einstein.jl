using FFTW
using FastTransforms

"""
    cheb2_quad_wts([TR=Float64], n::Integer)

Compute quadrature weights for Chebyshev points of the second kind (Clenshaw-Curtis quadrature).

# Arguments
- `TR`: Type parameter for the weights (e.g., Float64)
- `n`: Number of points

# Returns
- Vector of n weights for Clenshaw-Curtis quadrature

# Mathematical Background
For a function expressed in the Chebyshev basis:
```math
f(x) = \\sum_{k=0}^n c_k T_k(x)
```
The definite integral can be expressed as:
```math
\\int_{-1}^1 f(x)dx = \\mathbf{v}^T\\mathbf{c}
```
where ``\\mathbf{v}`` contains the integrals of Chebyshev polynomials:
```math
v_k = \\int_{-1}^1 T_k(x)dx = \\begin{cases}
\\frac{2}{1-k^2} & k \\text{ even} \\\\
0 & k \\text{ odd}
\\end{cases}
```

The quadrature weights ``\\mathbf{w}`` satisfy:
```math
\\int_{-1}^1 f(x)dx \\approx \\sum_{j=1}^n w_j f(x_j)
```

# Algorithm
Uses Waldvogel's algorithm (2006) with modifications by Nick Hale:
1. Compute exact integrals of even-indexed Chebyshev polynomials
2. Mirror the sequence for DCT via FFT
3. Apply inverse FFT
4. Adjust boundary weights

# Edge Cases
- n = 0: Returns empty vector
- n = 1: Returns [2.0]
- n ≥ 2: Returns full set of weights

# Examples
```julia
# Compute weights for 5-point quadrature
w = cheb2_quad_wts(5)

# Integrate sin(x) from -1 to 1
x = cheb2_grid(5)
f = sin.(x)
I = dot(w, f)  # ≈ 0
```

# References
1. [waldvogel2006fast](@citet*)
2. [Fast Clenshaw-Curtis Quadrature - File Exchange - MATLAB Central](https://www.mathworks.com/matlabcentral/fileexchange/6911-fast-clenshaw-curtis-quadrature)
"""
function cheb2_quad_wts(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    if n == 0
        return TR[]
    elseif n == 1
        return TR[2]
    end

    # Exact integrals of T_k (even k)
    c = [TR(2) / (1 - k^2) for k in [0; 2:2:(n - 1)]]

    # Mirror for DCT via FFT
    append!(c, c[floor(Int, n / 2):-1:2])

    # Compute weights via inverse FFT
    w = real.(ifft(c))

    # Adjust boundary weights
    w[1] /= 2
    append!(w, w[1])

    return w
end

function cheb2_quad_wts(n::TI) where {TI<:Integer}
    return cheb2_quad_wts(Float64, n)
end

export cheb2_quad_wts
