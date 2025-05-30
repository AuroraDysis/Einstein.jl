"""
Copyright (c) 2025 Zhen Zhong
Copyright (c) 2019-2020 Elias Jarlebring, Max Bennedich, Giampaolo Mele, Emil Ringh, Parikshit Upadhyaya

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
    A, E = qnm_pep_companion(TFC, pep::AbstractVector{<:AbstractMatrix}) where TFC <: Union{AbstractFloat, Complex{<:AbstractFloat}}

Linearizes a polynomial eigenvalue problem (PEP) and returns the companion form, as in the paper by Mehrmann and Voss.
More precisely, for a k-th degree PEP with n-by-n coefficient matrices,
this returns matrices A and E, both kn-by-kn, corresponding to the linearized problem
```math
Ax = λEx
```

# References
- [mehrmann2004nonlinear](@citet*)
- [NonlinearEigenproblems.jl/src/method_companion.jl at master · nep-pack/NonlinearEigenproblems.jl](https://github.com/nep-pack/NonlinearEigenproblems.jl/blob/master/src/method_companion.jl)
"""
function qnm_pep_companion(pep::AbstractVector{<:AbstractMatrix{TN}}) where {TN<:Number}
    @boundscheck begin
        @argcheck !isempty(pep) "pep must not be empty"
        @argcheck all(m -> size(m, 1) == size(pep[1], 1), pep) "all matrices must have the same size"
        @argcheck all(m -> size(m, 1) == size(m, 2), pep) "all matrices must be square"
    end

    is_sparse = issparse(pep[1])

    # Size of coefficient matrices
    n = size(pep[1], 1)

    # Degree of pep
    d = length(pep) - 1

    Iblock = kron(Eye{Int}(d - 1), Eye{Int}(n))

    if is_sparse
        A = spzeros(TN, n * d, n * d)
        E = spzeros(TN, d * n, d * n)
    else
        A = zeros(TN, n * d, n * d)
        E = zeros(TN, d * n, d * n)
    end

    # -- Construct A -- #

    # First row block of A
    for i in 1:d
        @.. A[1:n, ((i - 1) * n + 1):(i * n)] = -pep[d - i + 1]
    end

    # Lower part of A
    A[(n + 1):(d * n), 1:((d - 1) * n)] = Iblock

    # -- Construct E -- #

    # Fill block (1,1)
    E[1:n, 1:n] = pep[d + 1]

    #Fill all blocks on the diagonal with eye(n)
    E[(n + 1):(d * n), (n + 1):(d * n)] = Iblock

    return A, E
end

export qnm_pep_companion
