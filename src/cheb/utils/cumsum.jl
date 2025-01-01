"""
    cumsum(f::Vector{TR}) where {TR<:AbstractFloat}

Compute the indefinite integral of a function represented in Chebyshev basis.

# Mathematical Background
Given a Chebyshev series of length n:
```math
G(x) = \\sum_{r=0}^{n-1} c_r T_r(x)
```
its integral is represented with a series of length n+1:
```math
\\int G(x)\\,dx = \\sum_{r=0}^{n} b_r T_r(x)
```

# Arguments
- `f::Vector{TR}`: Coefficients ``c_r`` of the Chebyshev series

# Returns
- Vector of coefficients ``b_r`` for the integral, with length n+1

# Mathematical Details
The coefficients are computed as follows:
```math
\\begin{align*}
b_0 &= \\sum_{r=1}^{n} (-1)^{r+1} b_r \\quad \\text{(constant of integration)} \\\\
b_1 &= c_0 - \\frac{c_2}{2} \\\\
b_r &= \\frac{c_{r-1} - c_{r+1}}{2r} \\quad \\text{for } r > 1
\\end{align*}
```
where ``c_{n+1} = c_{n+2} = 0``

The constant term ``b_0`` is chosen to make ``f(-1) = 0``.

# Examples
```julia
using Test

# Suppose we have 3 Chebyshev coefficients:
f = [1.0, 2.0, 3.0]

# Perform the continuous sum:
f_new = cumsum(f)

@test length(f_new) == length(f) + 1
println("Original coefficients: ", f)
println("New coefficients: ", f_new)
```

# References
- Mason & Handscomb, "Chebyshev Polynomials", Chapman & Hall/CRC (2003), pp. 32-33.
- [chebfun/@chebtech/cumsum.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech/cumsum.m)

See also: [`diff`](@ref), [`sum`](@ref)
"""
function cumsum(f::Vector{TR}) where {TR<:AbstractFloat}
    # Copy input coefficients
    c = copy(f)
    n = length(c)

    # Pad with zeros (equivalent to c = [c; zeros(2, m)] in MATLAB with m=1)
    c = vcat(c, zeros(eltype(c), 2))

    # Prepare a new vector b (n+1 in length)
    b = zeros(eltype(c), n + 1)

    # Compute b[3] ... b[n+1] (1-based indexing in Julia):
    # b(3:n+1) = (c(2:n) - c(4:n+2)) ./ (2*(2:n)) in MATLAB
    for r in 2:n
        b[r + 1] = (c[r] - c[r + 2]) / (2 * r)
    end

    # b(2) = c(1) - c(3)/2
    b[2] = c[1] - c[3] / 2

    # Construct vector v where v = [1, -1, 1, -1, ...]
    v = ones(eltype(c), n)
    for i in 2:2:n
        v[i] = -1
    end

    # b(1) = v * b(2:end) (dots with appropriate indices)
    # This enforces the condition analogous to f(-1) = 0 from the MATLAB code.
    # If you want to further adjust the constant term for a specific boundary condition,
    # you can do so here by evaluating the polynomial at x=-1 and subtracting that value.
    b[1] = sum(v .* b[2:end])

    return b
end

@testset "cumsum" begin
    tol = 100 * eps()

    @testset "cheb1" begin
        n = 15
        f = cos.(cheb1_pts(n))
        f_coeffs = cheb1_vals2coeffs(f)
        If_coeffs = cumsum(f_coeffs)
        If = sin.(cheb1_pts(n)) .- sin(-1) # sin(x) - sin(-1) is the antiderivative of cos(x)
        If_coeffs_true = cheb1_vals2coeffs(If)
        @test norm(If_coeffs[1:end-1] .- If_coeffs_true, Inf) < tol
    end

    @testset "cheb2" begin
        n = 15
        f = cos.(cheb2_pts(n))
        f_coeffs = cheb2_vals2coeffs(f)
        If_coeffs = cumsum(f_coeffs)
        If = sin.(cheb2_pts(n)) .- sin(-1) # sin(x) - sin(-1) is the antiderivative of cos(x)
        If_coeffs_true = cheb2_vals2coeffs(If)
        @test norm(If_coeffs[1:end-1] .- If_coeffs_true, Inf) < tol
    end
end
