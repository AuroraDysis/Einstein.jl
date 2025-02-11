"""
    cheb1_amat([TF=Float64], n::Integer) where {TF<:AbstractFloat}

Construct the analysis matrix A that transforms function values at Chebyshev points of the 1st kind to Chebyshev coefficients.

# Arguments
- `TF`: Element type (defaults to Float64)
- `n`: Number of points/coefficients
"""
function cheb1_amat(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    A = Array{TF,2}(undef, n, n)
    op = Cheb1Vals2CoeffsOp{TF}(n)
    @inbounds for i in 1:n
        A[:, i] = op(OneElement(one(TF), i, n))
    end
    return A
end

function cheb1_amat(n::Integer)
    return cheb1_amat(Float64, n)
end

"""
    cheb1_smat([TF=Float64], n::Integer) where {TF<:AbstractFloat}

Construct the synthesis matrix S that transforms Chebyshev coefficients to function values at Chebyshev points of the 1st kind.

# Arguments
- `TF`: Element type (defaults to Float64)
- `n`: Number of points/coefficients
"""
function cheb1_smat(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    S = Array{TF,2}(undef, n, n)
    op = Cheb1Coeffs2ValsOp{TF}(n)
    @inbounds for i in 1:n
        S[:, i] = op(OneElement(one(TF), i, n))
    end
    return S
end

function cheb1_smat(n::Integer)
    return cheb1_smat(Float64, n)
end

export cheb1_amat, cheb1_smat
