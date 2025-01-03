"""
    cheb1_amat([TR=Float64], n::Integer)

Construct the analysis matrix A that transforms function values at Chebyshev points of the 1st kind to Chebyshev coefficients.

# Arguments
- `TR`: Element type (defaults to Float64)
- `n`: Number of points/coefficients
"""
function cheb1_amat(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    A = Array{TR,2}(undef, n, n)
    op = Cheb1Vals2CoeffsOp(TR, n)
    @inbounds for i in 1:n
        A[:, i] = op(OneElement(one(TR), i, n))
    end
    return A
end

function cheb1_amat(n::TI) where {TI<:Integer}
    return cheb1_amat(Float64, n)
end

"""
    cheb1_smat([TR=Float64], n::Integer)

Construct the synthesis matrix S that transforms Chebyshev coefficients to function values at Chebyshev points of the 1st kind.

# Arguments
- `TR`: Element type (defaults to Float64)
- `n`: Number of points/coefficients
"""
function cheb1_smat(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    S = Array{TR,2}(undef, n, n)
    op = Cheb1Coeffs2ValsOp(TR, n)
    @inbounds for i in 1:n
        S[:, i] = op(OneElement(one(TR), i, n))
    end
    return S
end

function cheb1_smat(n::TI) where {TI<:Integer}
    return cheb1_smat(Float64, n)
end

export cheb1_amat, cheb1_smat
