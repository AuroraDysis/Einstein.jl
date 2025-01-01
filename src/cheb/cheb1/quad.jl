using FFTW
using FastTransforms

"""
    cheb1_quad_wts([TR=Float64], n::Integer)

Compute quadrature weights for Chebyshev points of the first kind.

# Arguments
- `TR`: Type parameter for the weights (e.g., Float64)
- `n`: Number of points

# Returns
- Vector of n weights for Chebyshev quadrature

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
w = cheb1_quad_wts(5)

# Integrate sin(x) from -1 to 1
x = cheb1_pts(5)
f = sin.(x)
I = dot(w, f)  # ≈ 0
```

# References
1. [waldvogel2006fast](@citet*)
2. [Fast Clenshaw-Curtis Quadrature - File Exchange - MATLAB Central](https://www.mathworks.com/matlabcentral/fileexchange/6911-fast-clenshaw-curtis-quadrature)
"""
function cheb1_quad_wts(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    # Handle the special cases:
    if n == 0
        return TR[]
    elseif n == 1
        return TR[2]
    end

    # Preallocate the array m for the moments
    evens = 2:2:(n - 1)
    m = zeros(TR, length(evens) + 1)
    @inbounds m[1] = 2 / one(TR)  # Corresponds to k=0
    @inbounds for (i, k) in enumerate(evens)
        # m[i+1] = 1 - k^2
        m[i + 1] = 2 / (one(TR) - k^2)
    end

    # Preallocate the coefficient array for FFT
    c = zeros(Complex{TR}, n)
    c[1:length(m)] .= m

    # Fill the coefficient array based on parity of n
    if isodd(n)
        # For odd n: c = [m..., -reverse(m[(n+1)/2 : 2])]
        c[(length(m) + 1):end] .= -@view(m[div(n + 1, 2):-1:2])
    else
        # For even n: c = [m..., 0, -reverse(m[n/2 : 2])]
        c[(length(m) + 2):end] .= -@view(m[div(n, 2):-1:2])
    end

    # Multiply by rotation factors exp(1im*(0:n-1)*π/n)
    im_pi_over_n = im * convert(TR, π) / n
    @inbounds for k in 0:(n - 1)
        c[k + 1] *= exp(k * im_pi_over_n)
    end

    # Compute inverse FFT in-place and take the real part for the weights
    ifft!(c)
    w = real.(c)

    return w
end

function cheb1_quad_wts(n::TI) where {TI<:Integer}
    return cheb1_quad_wts(Float64, n)
end

export cheb1_quad_wts

@testset "cheb1_quad_wts" begin
    # Test n=0 case
    @test cheb1_quad_wts(0) == Float64[]

    # Test n=1 case
    @test cheb1_quad_wts(1) ≈ [2.0]

    # Test n=5 case
    w5 = cheb1_quad_wts(5)
    @test w5 ≈ [
        0.167781228466683,
        0.525552104866650,
        0.613333333333333,
        0.525552104866650,
        0.167781228466684,
    ]

    w6 = cheb1_quad_wts(6)
    @test w6 ≈ [
        0.118661021381236,
        0.377777777777778,
        0.503561200840986,
        0.503561200840986,
        0.377777777777778,
        0.118661021381236,
    ]
end
