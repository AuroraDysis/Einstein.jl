"""
    cheb2_amat([T=Float64], n::Integer) where {T<:AbstractFloat}

Construct the analysis matrix A that transforms function values at Chebyshev points of the 2nd kind to Chebyshev coefficients.

# Arguments
- `T`: Element type (defaults to Float64)
- `n`: Number of points/coefficients
"""
function cheb2_amat(::Type{T}, n::Integer) where {T<:AbstractFloat}
    A = Array{T,2}(undef, n, n)
    op = Cheb2Vals2CoeffsOp{T}(n)
    @inbounds for i in 1:n
        A[:, i] = op(OneElement(one(T), i, n))
    end
    return A
end

function cheb2_amat(n::Integer)
    return cheb2_amat(Float64, n)
end

"""
    cheb1_smat([T=Float64], n::Integer) where {T<:AbstractFloat}

Construct the synthesis matrix S that transforms Chebyshev coefficients to function values at Chebyshev points of the 2nd kind.

# Arguments
- `T`: Element type (defaults to Float64)
- `n`: Number of points/coefficients
"""
function cheb2_smat(::Type{T}, n::Integer) where {T<:AbstractFloat}
    S = Array{T,2}(undef, n, n)
    op = Cheb2Coeffs2ValsOp{T}(n)
    @inbounds for i in 1:n
        S[:, i] = op(OneElement(one(T), i, n))
    end
    return S
end

function cheb2_smat(n::Integer)
    return cheb2_smat(Float64, n)
end

export cheb2_amat, cheb2_smat
