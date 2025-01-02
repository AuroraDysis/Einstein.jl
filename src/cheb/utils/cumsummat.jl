
"""
    cheb_coeffs_cumsummat([TR=Float64], n::TI) where {TR<:AbstractFloat,TI<:Integer}
    cheb_coeffs_cumsummat([TR=Float64], n::TI, x_min::TR, x_max::TR) where {TR<:AbstractFloat,TI<:Integer}

Generate the Chebyshev coefficient integration matrix
that maps Chebyshev coefficients to the coefficients
of the integral of the interpolating polynomial.

# Arguments
- `TR`: Type parameter for the matrix elements (e.g., Float64)
- `n`: Size of the matrix (n×n)
- `x_min`: (Optional) Lower bound of the integration interval
- `x_max`: (Optional) Upper bound of the integration interval

# References
- [chebfun/@chebcolloc1/chebcolloc1.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebcolloc1/chebcolloc1.m)
- [chebfun/@chebcolloc2/chebcolloc2.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebcolloc2/chebcolloc2.m)
"""
function cheb_coeffs_cumsummat(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
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

function cheb_coeffs_cumsummat(n::TI) where {TI<:Integer}
    return cheb_coeffs_cumsummat(Float64, n)
end

# Second method documentation is inherited from the main docstring
function cheb_coeffs_cumsummat(
    ::Type{TR}, n::TI, x_min::TR, x_max::TR
) where {TR<:AbstractFloat,TI<:Integer}
    B = cheb_coeffs_cumsummat(TR, n)
    B .*= (x_max - x_min) / 2
    return B
end

function cheb_coeffs_cumsummat(n::TI, x_min::Float64, x_max::Float64) where {TI<:Integer}
    return cheb_coeffs_cumsummat(Float64, n, x_min, x_max)
end
