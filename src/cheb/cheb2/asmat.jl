"""
    cheb2_amat([TF=Float64], n::Integer) where {TF<:AbstractFloat}

Construct the analysis matrix A that transforms function values at Chebyshev points of the 2nd kind to Chebyshev coefficients.

# Arguments
- `TF`: Element type (defaults to Float64)
- `n`: Number of points/coefficients
"""
function cheb2_amat(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    A = Array{TF,2}(undef, n, n)
    op = Cheb2Vals2CoeffsOp{TF}(n)
    @inbounds for i in 1:n
        A[:, i] = op(OneElement(one(TF), i, n))
    end
    return A
end

function cheb2_amat(n::Integer)
    return cheb2_amat(Float64, n)
end

"""
    cheb2_smat([TF=Float64], n::Integer) where {TF<:AbstractFloat}

Construct the synthesis matrix S that transforms Chebyshev coefficients to function values at Chebyshev points of the 2nd kind.

# Arguments
- `TF`: Element type (defaults to Float64)
- `n`: Number of points/coefficients
"""
function cheb2_smat(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    S = Array{TF,2}(undef, n, n)
    op = ChebyshevSecondKindSynthesis{TF}(n)
    @inbounds for i in 1:n
        S[:, i] = op(OneElement(one(TF), i, n))
    end
    return S
end

function cheb2_smat(n::Integer)
    return cheb2_smat(Float64, n)
end

export cheb2_amat, cheb2_smat
