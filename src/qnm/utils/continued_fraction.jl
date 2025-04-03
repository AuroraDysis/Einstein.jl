"""
MIT License

Copyright (c) 2025 Zhen Zhong
Copyright (c) 2019 Leo C. Stein

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

@doc raw"""
    continued_fraction_lentz([T=Float64], a::Function, b::Function, tol::T, min_iter::Integer, max_iter::Integer) where {T<:AbstractFloat}

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
continued_fraction_lentz(Float64, a, b, 10*eps(Float64), 50, 1000)
```

## Compute Golden Ratio

```math
\phi = 1 + \frac{1}{1 + \frac{1}{1 + \frac{1}{1 + \cdots}}} \approx 1.618033988749895
```

```julia
a(i) = 1
b(i) = 1
continued_fraction_lentz(Float64, a, b, 10*eps(Float64), 50, 1000)
```

# References
- [qnm/qnm/contfrac.py at master · duetosymmetry/qnm](https://github.com/duetosymmetry/qnm/blob/master/qnm/contfrac.py)
- [press2007numerical, Stein:2019mop, lentz1976generating, thompson1986coulomb, DLMF_3103_online](@cite)
"""
function continued_fraction_lentz(
    ::Type{T}, a::Function, b::Function, tol::T, min_iter::Integer, max_iter::Integer
) where {T<:AbstractFloat}
    tiny = eps(T)^2

    fᵢ₋₁ = b(0)
    if iszero(fᵢ₋₁)
        fᵢ₋₁ = tiny
    end
    Cᵢ₋₁ = fᵢ₋₁
    Dᵢ₋₁ = zero(fᵢ₋₁)

    fᵢ = fᵢ₋₁ # use fᵢ to store the result
    errorᵢ = typemax(T) # use errorᵢ to store the error

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
                "continued_fraction_lentz: Continued fraction did not converge after $max_iter iterations",
            ),
        )
    end

    return fᵢ, errorᵢ, i
end

function continued_fraction_lentz(
    a::Function, b::Function, tol::Float64, min_iter::Integer, max_iter::Integer
)
    return continued_fraction_lentz(Float64, a, b, tol, min_iter, max_iter)
end

export continued_fraction_lentz
