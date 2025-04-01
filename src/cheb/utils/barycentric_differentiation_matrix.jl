"""
    barycentric_differentiation_matrix(x, w, k, t)

Compute the barycentric differentiation matrix.

# Arguments
- `x::AbstractVector{TF}` : Vector of interpolation points
- `w::AbstractVector{TF}` : Barycentric weights of the interpolation points
- `k::Integer` : Order of the derivative (default: 1)
- `t::AbstractVector{TF}` : Vector of angles (default: empty)

# References:
- [chebfun/@chebcolloc/baryDiffMat.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebcolloc/baryDiffMat.m)
"""
function barycentric_differentiation_matrix(
    x::AbstractVector{TF}, w::AbstractVector{TF}, k::Integer=1, t::AbstractVector{TF}=TF[]
) where {TF<:AbstractFloat}
    @boundscheck begin
        @argcheck length(x) == length(w) "x and w must have the same length"
    end

    n = length(x)

    # Handle trivial cases:
    if n == 0
        return Matrix{TF}(undef, 0, 0)
    elseif n == 1
        # Derivative of a constant is 0
        return zeros(TF, 1, 1)
    end

    # 0-th derivative is the identity matrix
    if k == 0
        return Matrix{TF}(I, n, n)
    end

    # Prepare pairwise difference matrix Dx and pairwise weight ratio Dw:
    # We'll construct them so that Dx[i,j] = x[i] - x[j], Dw[i,j] = w[j]/w[i].
    Dx = zeros(TF, n, n)
    Dw = zeros(TF, n, n)

    # If t is provided, use the sine-based formula for differences (see reference [4]):
    if length(t) > 0
        # In MATLAB code, t was optionally passed in a reversed order. 
        # Here we assume user passes the correct array. 
        # We'll mimic the "flipud" approach if needed.
        t = collect(reverse(t))
        # Build matrix using: 2*sin((t[i]+t[j])/2)*sin((t[i]-t[j])/2)
        half = one(TF) / 2
        for i in 1:n
            for j in 1:n
                if i == j
                    Dx[i, j] = 1  # set the diagonal to 1 (placeholder)
                else
                    Dx[i, j] = 2 * sin((t[i] + t[j]) * half) * sin((t[i] - t[j]) * half)
                end
            end
        end

        # Instead, rotate Dx by 180 degrees, then flip signs above diagonal:
        DxRot = reverse(reverse(Dx; dims=1); dims=2)
        for i in 1:n
            for j in 1:n
                if j > i
                    Dx[i, j] = -DxRot[i, j]
                end
            end
        end

        # Finally, set the diagonal to 1:
        for i in 1:n
            Dx[i, i] = 1
        end
    else
        # Standard difference (x[i] - x[j]):
        for i in 1:n
            for j in 1:n
                if i == j
                    Dx[i, j] = 1
                else
                    Dx[i, j] = x[i] - x[j]
                end
            end
        end
    end

    # Build Dw = w[j] / w[i], with zero on diagonal:
    for i in 1:n
        for j in 1:n
            if i == j
                Dw[i, j] = 0
            else
                Dw[i, j] = w[j] / w[i]
            end
        end
    end

    # Compute reciprocal of Dx:
    Dxi = one(TF) ./ Dx  # elementwise reciprocal
    # Construct initial derivative matrix:
    D = Dw .* Dxi

    # Negative sum trick on the diagonal: D(i,i) = -sum of row i
    # This ensures each row sums to zero.
    negative_sum_trick!(D)

    # Enforce a symmetry fix for even n on the diagonal entries:
    # (Matches the "Forcing symmetry for even n" from the original code.)
    # The original code effectively flips the sign for half of the diagonal indices for even n.
    # For simplicity, replicate that behavior exactly:
    halfN = div(n, 2)
    # Indices from end down to n-floor(n/2)+1:
    # That is from n down to n-halfN+1
    Ddiag = collect(diag(D))
    Ddiag[end:-1:(n - halfN + 1)] = -Ddiag[1:halfN]
    for i in 1:n
        D[i, i] = Ddiag[i]
    end

    # If k=1, we're done:
    if k == 1
        return D
    end

    # For k=2:
    # D <- 2 * D .* (repeat(diag(D), 1, n) - Dxi)
    # and then negative sum trick
    D2 = copy(D)
    for i in 1:n
        # diag(D) for row i is D[i,i]
        dii = D[i, i]
        for j in 1:n
            D2[i, j] = 2 * D[i, j] * (dii - Dxi[i, j])
        end
    end
    # negative sum trick
    negative_sum_trick!(D2)

    if k == 2
        return D2
    end

    # For k >= 3:
    # We'll iteratively build D for nth derivative:
    D_current = copy(D2)
    for n in 3:k
        # D <- n * Dxi .* (Dw .* repeat(diag(D_current), 1, n) - D_current)
        D_next = similar(D_current)
        for i in 1:n
            dii = D_current[i, i]
            for j in 1:n
                D_next[i, j] = n * Dxi[i, j] * (Dw[i, j] * dii - D_current[i, j])
            end
        end
        # negative sum trick
        negative_sum_trick!(D_next)
        D_current = D_next
    end

    return D_current
end

function negative_sum_trick!(D::Matrix{TF}) where {TF<:AbstractFloat}
    D[diagind(D)] .= zero(TF)
    return D[diagind(D)] .-= sum(D; dims=2)
end

export barycentric_differentiation_matrix
