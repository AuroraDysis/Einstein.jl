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
    qnm_schw_expansion_dolan_ottewill(::Type{TF}, s::TI, l::TI, n::TI) where {TF<:AbstractFloat,TI<:Integer}

Compute the high ``\ell`` asymptotic expansion of Schwarzschild QNM frequency, using the method of Dolan and Ottewill [Dolan:2009nk](@cite).

The QNM frequency can be written as an expansion in inverse powers of ``L=\ell+\frac{1}{2}``
```math
\omega_{l n}=\varpi_{-1}^{(n)} L+\varpi_0^{(n)}+\varpi_1^{(n)} L^{-1}+\varpi_2^{(n)} L^{-2}+\ldots
```

The lowest expansion coefficients for arbitrary spin ``\beta=1-s^2`` and arbitrary overtone number ``n`` are given by
```math
\begin{aligned}
& \sqrt{27} \varpi_{-1}^{(n)}=1 \\
& \sqrt{27} \varpi_0^{(n)}=-i N \\
& \sqrt{27} \varpi_1^{(n)}=\frac{\beta}{3}-\frac{5 N^2}{36}-\frac{115}{432} \\
& \sqrt{27} \varpi_2^{(n)}=-i N\left[\frac{\beta}{9}+\frac{235 N^2}{3888}-\frac{1415}{15552}\right] \\
& \sqrt{27} \varpi_3^{(n)}=-\frac{\beta^2}{27}+\frac{204 N^2+211}{3888} \beta+\frac{854160 N^4-1664760 N^2-776939}{40310784} \\
& \sqrt{27} \varpi_4^{(n)}=i N\left[\frac{\beta^2}{27}+\frac{1100 N^2-2719}{46656} \beta+\frac{11273136 N^4-52753800 N^2+66480535}{2902376448}\right]
\end{aligned}
```

# Arguments

- `TF`: Type of the floating-point number.
- `s`: Spin weight of the field of interest.
- `l`: Multipole number of interest.
- `n`: Overtone number of interest.

# References
- [Dolan:2009nk](@cite)
"""
function qnm_schw_expansion_dolan_ottewill(
    ::Type{TF}, s::Integer, l::Integer, n::Integer
) where {TF<:AbstractFloat}
    L = l + one(TF) / 2
    N = n + one(TF) / 2
    β = one(TF) - s * s

    ϖ₋₁ = one(TF)
    ϖ₀ = -1im * N
    ϖ₁ = β / 3 - 5 * N * N / 36 - TF(115) / 432
    ϖ₂ = -1im * N * (β / 9 + 235 * N * N / 3888 - TF(1415) / 15552)
    ϖ₃ =
        -β * β / 27 +
        (204 * N * N + 211) / 3888 * β +
        (854160 * N^4 - 1664760 * N * N - 776939) / 40310784
    ϖ₄ =
        -1im *
        N *
        (
            β * β / 27 +
            (1100 * N * N - 2719) / 46656 * β +
            (11273136 * N^4 - 52753800 * N * N + 66480535) / 2902376448
        )
    ϖᵢ = SA[ϖ₋₁, ϖ₀, ϖ₁, ϖ₂, ϖ₃, ϖ₄]
    ω = zero(Complex{TF})
    @inbounds for n in -1:4
        ω += ϖᵢ[n + 2] * L^(-n)
    end

    return ω / sqrt(TF(27))
end

export qnm_schw_expansion_dolan_ottewill
