@doc raw"""
    contfrac_lentz([TR=Float64], a::TF1, b::TF2, tol::TR, min_iter::TI, max_iter::TI) where { T <: AbstractFloat, TF1 <: Function, TF2 <: Function }

Compute the continued fraction
```math
f(x)=b_0+\frac{a_1}{b_1+\frac{a_2}{b_2+\frac{a_3}{b_3+\frac{a_4}{b_4+\frac{a_5}{b_5+\cdots}}}}}
```
using modified Lentz's method. Translated from [duetosymmetry/qnm](https://github.com/duetosymmetry/qnm).

# Arguments
- `a`: A function that returns the `aᵢ` terms.
- `b`: A function that returns the `bᵢ` terms.
- `tol`: The tolerance for convergence.
- `min_iter`: The minimum number of iterations to perform.
- `max_iter`: The maximum number of iterations to perform.

# Returns
- `fᵢ`: The value of the continued fraction.
- `errorᵢ`: The estimated error.
- `i`: The number of iterations performed.

# Examples

## Compute the square root of two using continued fractions

```math
\sqrt{2} = 1 + \frac{1}{2 + \frac{1}{2 + \frac{1}{2 + \frac{1}{2 + \cdots}}}} \approx 1.414213562373095
```

```julia
a(i) = 1
b(i) = i == 0 ? 1 : 2
contfrac_lentz(Float64, a, b, 10*eps(Float64), 50, 1000)
```

## Compute Golden Ratio

```math
\phi = 1 + \frac{1}{1 + \frac{1}{1 + \frac{1}{1 + \cdots}}} \approx 1.618033988749895
```

```julia
a(i) = 1
b(i) = 1
contfrac_lentz(Float64, a, b, 10*eps(Float64), 50, 1000)
```

# References
- [qnm/qnm/contfrac.py at master · duetosymmetry/qnm](https://github.com/duetosymmetry/qnm/blob/master/qnm/contfrac.py)
- [press2007numerical, Stein:2019mop, lentz1976generating, thompson1986coulomb, DLMF§3103:online](@cite)
"""
function contfrac_lentz(
    ::Type{TR}, a::TF1, b::TF2, tol::TR, min_iter::TI, max_iter::TI
) where {TR<:AbstractFloat,TF1<:Function,TF2<:Function,TI<:Integer}
    tiny = eps(TR)^2

    fᵢ₋₁ = b(0)
    if iszero(fᵢ₋₁)
        fᵢ₋₁ = tiny
    end
    Cᵢ₋₁ = fᵢ₋₁
    Dᵢ₋₁ = zero(fᵢ₋₁)

    fᵢ = fᵢ₋₁ # use fᵢ to store the result
    errorᵢ = typemax(TR) # use errorᵢ to store the error

    i = 1
    converged = false
    while i <= max_iter
        aᵢ = a(i)
        bᵢ = b(i)

        Dᵢ = bᵢ + aᵢ * Dᵢ₋₁
        if iszero(Dᵢ)
            Dᵢ = tiny
        end

        Cᵢ = bᵢ + aᵢ / Cᵢ₋₁
        if iszero(Cᵢ)
            Cᵢ = tiny
        end

        Dᵢ = inv(Dᵢ)
        Δᵢ = Cᵢ * Dᵢ
        fᵢ = fᵢ₋₁ * Δᵢ
        errorᵢ = abs(Δᵢ - 1)

        if i >= min_iter && errorᵢ < tol
            converged = true
            break
        end

        i += 1
        Dᵢ₋₁ = Dᵢ
        Cᵢ₋₁ = Cᵢ
        fᵢ₋₁ = fᵢ
    end

    if !converged
        throw(
            ConvergenceError(
                "contfrac_lentz: Continued fraction did not converge after $max_iter iterations",
            ),
        )
    end

    return fᵢ, errorᵢ, i
end

export contfrac_lentz
