"""
    ultra_diffmat(n::TI, m::TI) where {TI<:Integer}

Differentiation matrices for ultraspherical spectral method that takes N Chebyshev
coefficients and returns \$n C^{(m)}\$ coefficients that represent the derivative
of the Chebyshev series. Here, \$C^{(k)}\$ is the ultraspherical polynomial basis
with parameter \$k\$.

# Arguments
- `n::TI`: Number of points
- `m::TI`: Order of differentiation

# References
- [chebfun/@ultraS/diffmat.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40ultraS/diffmat.m)
"""
function ultra_diffmat(n::TI, m::TI) where {TI<:Integer}
    @argcheck n >= 1 "n must be positive"
    @argcheck m >= 0 "m must be nonnegative"

    nm1 = n - 1
    if m > 0
        D = spdiagm(1 => 1:nm1)
        for s in 1:(m - 1)
            D = spdiagm(1 => 2s * ones(TI, nm1)) * D
        end
    else
        D = sparse(TI, I, nm1, nm1)
    end

    return D
end

export ultra_diffmat

@testset "ultra_diffmat" begin
    for n in 2:10
        for m in 1:3
            @test isapprox(ultra_diffmat(n, m), Derivative(Chebyshev(), m)[1:n, 1:n])
        end
    end
end
