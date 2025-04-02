"""
    barycentric_differentiation_matrix(x::AbstractVector{TF}, w::AbstractVector{TF}, k::Integer=1, t::Union{AbstractVector{TF},Nothing}=nothing) where {TF<:AbstractFloat}

Compute the barycentric differentiation matrix.

# Arguments
- `x::AbstractVector{TF}` : Vector of interpolation points
- `w::AbstractVector{TF}` : Barycentric weights of the interpolation points
- `k::Integer` : Order of the derivative (default: 1)
- `t::AbstractVector{TF}` : Vector of angles (default: nothing)

# References:
- [chebfun/@chebcolloc/baryDiffMat.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebcolloc/baryDiffMat.m)
- [Baltensperger2003](@citet*)
"""
function barycentric_differentiation_matrix(
    x::AbstractVector{TF},
    w::AbstractVector{TF},
    k::Integer=1,
    t::Union{AbstractVector{TF},Nothing}=nothing,
) where {TF<:AbstractFloat}
    @argcheck length(x) == length(w) "x and w must have the same length"

    n = length(x)

    # handle trivial cases:
    if n == 0
        return Matrix{TF}(undef, 0, 0)
    elseif n == 1
        # derivative of a constant is 0
        return zeros(TF, 1, 1)
    end

    # 0-th derivative is the identity matrix
    if k == 0
        return Matrix{TF}(I, n, n)
    end

    # Dx[i,j] = x[i] - x[j]
    Dx = zeros(TF, n, n)
    diag_idx = diagind(Dx)

    if !isnothing(t)
        @argcheck length(t) == n "t must have the same length as x"

        # use trigonometric identities as described in Baltensperger2003
        tr = @view t[end:-1:1]
        for j in 1:n, i in 1:(n - j + 1)
            Dx[i, j] = 2 * sin((tr[i] + tr[j]) / 2) * sin((tr[i] - tr[j]) / 2)
        end
    else
        # standard pairwise differences: Dx[i, j] = x[i] - x[j]
        for j in 1:n, i in 1:(n - j + 1)
            Dx[i, j] = x[i] - x[j]
        end
    end

    # flipping trick, D_{N-k, N-j}=-D_{k j}
    for j in 1:n, i in (n - j + 2):n
        Dx[i, j] = -Dx[n - i + 1, n - j + 1]
    end

    # build Dw = w[j] / w[i], with zero on diagonal:
    Dw = zeros(TF, n, n)
    for i in 1:n, j in 1:n
        Dw[i, j] = w[j] / w[i]
    end
    Dw[diag_idx] .= zero(TF)

    # k = 1
    Dx[diag_idx] .= one(TF) # Set diagonal to 1 temporarily to avoid division by zero
    Dxi = one(TF) ./ Dx
    D = Dw .* Dxi

    negative_sum_trick!(D)

    # enforce a symmetry fix for even n on the diagonal entries:
    half_n = div(n, 2)
    D[diag_idx[end:-1:(n - half_n + 1)]] .= -D[diag_idx[1:half_n]]

    if k == 1
        return D
    end

    # k = 2
    for i in 1:n
        Dii = D[i, i]
        for j in 1:n
            D[i, j] = 2 * D[i, j] * (Dii - Dxi[i, j])
        end
    end

    negative_sum_trick!(D)

    if k == 2
        return D
    end

    # For k >= 3
    for n in 3:k
        for i in 1:n
            Dii = D[i, i]
            for j in 1:n
                D[i, j] = n * Dxi[i, j] * (Dw[i, j] * Dii - D[i, j])
            end
        end

        negative_sum_trick!(D)
    end

    return D
end

# D_{k k}=-\sum_{\substack{j=0 \\ j \neq k}}^N D_{k j}
function negative_sum_trick!(D::Matrix{TF}) where {TF<:AbstractFloat}
    n = size(D, 1)

    for i in 1:n
        Dii = zero(TF)
        for j in 1:(i - 1)
            Dii += D[i, j]
        end
        for j in (i + 1):n
            Dii += D[i, j]
        end
        D[i, i] = -Dii
    end
end

export barycentric_differentiation_matrix
