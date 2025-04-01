"""
The MIT License (MIT)

Copyright (c) 2025 Zhen Zhong
Copyright (c) 2012 Jeff Bezanson, Jeffrey Sarnoff, and other contributors

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

"""
    dot_xsum(x::StridedVector{T}, y::StridedVector{T}) where {T<:Real}

Compute the dot product of two vectors using extended precision accumulation.
Uses the `xsum` package for improved numerical accuracy.

# Note
- Very fast for large vectors, but a bit slower than Kahan summation for small vectors.
- Both input vectors must have the same length, which is not checked for performance reasons.

# References
- [neal2015fastexactsummationusing](@citet*)
- [JuliaMath/Xsum.jl](https://github.com/JuliaMath/Xsum.jl)
- [Radford Neal / xsum · GitLab](https://gitlab.com/radfordneal/xsum)
"""
function dot_xsum(x::StridedVector{Float64}, y::StridedVector{Float64})
    acc = XAccumulator(zero(Float64))
    @inbounds for i in eachindex(x)
        accumulate!(acc, x[i] * y[i])
    end
    return float(acc)
end

"""
    dot_kahan(v1::StridedVector{T}, v2::StridedVector{T}) where {T<:Number}

Neumaier's variant of Kahan summation algorithm to reduce numerical errors.

# Note
- Slower than `dot_xsum` for large vectors, but faster for small vectors.
- Similar performance to `dot_kahan`
- Uses loop unrolling for better performance while maintaining Kahan summation's
numerical stability. Processes two elements per iteration when possible.

# References
- [JuliaMath/KahanSummation.jl](https://github.com/JuliaMath/KahanSummation.jl)
"""
function dot_kahan(v1::StridedVector{T}, v2::StridedVector{T}) where {T<:Number}
    @inbounds begin
        n = length(v1)
        c = zero(T)
        s = v1[1] * v2[1] - c
        for i in 2:n
            si = v1[i] * v2[i]
            t = s + si
            if abs(s) >= abs(si)
                c -= ((s - t) + si)
            else
                c -= ((si - t) + s)
            end
            s = t
        end
        return s - c
    end
end

export dot_xsum, dot_kahan
