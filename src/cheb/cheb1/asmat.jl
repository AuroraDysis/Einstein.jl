"""
    cheb1_amat([T=Float64], n::Integer) where {T<:AbstractFloat}

Construct the analysis matrix A that transforms function values at Chebyshev points of the 1st kind to Chebyshev coefficients.

# Arguments
- `T`: Element type (defaults to Float64)
- `n`: Number of points/coefficients
"""
function cheb1_amat(::Type{T}, n::Integer) where {T<:AbstractFloat}
    A = Array{T,2}(undef, n, n)
    op = Cheb1Vals2CoeffsOp{T}(n)
    @inbounds for i in 1:n
        A[:, i] = op(OneElement(one(T), i, n))
    end
    return A
end

function cheb1_amat(n::Integer)
    return cheb1_amat(Float64, n)
end

"""
    cheb1_smat([T=Float64], n::Integer)

Construct the synthesis matrix S that transforms Chebyshev coefficients to function values at Chebyshev points of the 1st kind.

# Arguments
- `T`: Element type (defaults to Float64)
- `n`: Number of points/coefficients
"""
function cheb1_smat(::Type{T}, n::Integer) where {T<:AbstractFloat}
    S = Array{T,2}(undef, n, n)
    op = Cheb1Coeffs2ValsOp{T}(n)
    @inbounds for i in 1:n
        S[:, i] = op(OneElement(one(T), i, n))
    end
    return S
end

function cheb1_smat(n::Integer)
    return cheb1_smat(Float64, n)
end

export cheb1_amat, cheb1_smat
