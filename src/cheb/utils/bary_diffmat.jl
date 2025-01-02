"""
    bary_diffmat(x; w=nothing, k=1, t=nothing)

Compute the barycentric differentiation matrix for points in the vector `x`.

• `x` is a vector of distinct interpolation points.
• `w` is an vector of barycentric weights associated with `x`. 
• `k` is the order of the derivative to be approximated (default is 1).
• `t` is an optional vector (e.g., `acos.(x)`) that may improve numerical stability for
  certain sets of points (see references [4]).

Returns the matrix `D` such that `D * f_values ≈ f'` (or higher derivatives) at the points in `x`,
based on the barycentric interpolation formula.

References:
[1] Schneider, C., & Werner, W. (1986). Some new aspects of rational interpolation. Mathematics of Computation, 47(176), 285–299.
[2] Welfert, B. D. (1997). Generation of pseudospectral matrices I. SIAM Journal on Numerical Analysis, 34(4), 1640–1657.
[3] Tee, T. W. (2006). An adaptive rational spectral method for differential equations with rapidly varying solutions (DPhil Thesis). University of Oxford.
[4] Baltensperger, R., & Trummer, M. R. (2003). Spectral Differencing with a Twist. SIAM Journal on Scientific Computing, 24(5), 1465–1487.
[5] [chebfun/@chebcolloc/baryDiffMat.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebcolloc/baryDiffMat.m)
"""
function bary_diffmat(
    x::VT1, w::VT2, k::TI=1, t::Union{VT3,Nothing}=nothing
) where {
    TR<:AbstractFloat,
    VT1<:AbstractVector{TR},
    VT2<:AbstractVector{TR},
    VT3<:AbstractVector{TR},
    TI<:Integer,
}
    n = length(x)

    # Handle trivial cases:
    if n == 0
        return Matrix{TR}(undef, 0, 0)
    elseif n == 1
        return [zero(TR)]  # Single point -> derivative matrix is [0]
    end

    # Check for length agreement:
    if length(w) != n
        error("Length of x must match the length of w.")
    end

    # If k=0, return identity:
    if k == 0
        return Matrix{TR}(I, n, n)
    end

    # Prepare pairwise difference matrix Dx and pairwise weight ratio Dw:
    # We'll construct them so that Dx[i,j] = x[i] - x[j], Dw[i,j] = w[j]/w[i].
    Dx = zeros(TR, n, n)
    Dw = zeros(TR, n, n)

    # If t is provided, use the sine-based formula for differences (see reference [4]):
    if !isnothing(t)
        # In MATLAB code, t was optionally passed in a reversed order. 
        # Here we assume user passes the correct array. 
        # We'll mimic the "flipud" approach if needed.
        t = collect(reverse(t))
        # Build matrix using: 2*sin((t[i]+t[j])/2)*sin((t[i]-t[j])/2)
        half = one(TR) / 2
        for i in 1:n
            for j in 1:n
                if i == j
                    Dx[i, j] = 1  # set the diagonal to 1 (placeholder)
                else
                    Dx[i, j] = 2 * sin((t[i] + t[j]) * half) * sin((t[i] - t[j]) * half)
                end
            end
        end

        # Implement the flipping trick to enforce symmetry:
        # We'll rotate Dx 180 degrees and flip signs for upper/lower triangles.
        for i in 1:n
            for j in 1:n
                if j > i
                    # swap with bottom-left
                    Dx[i, j] = -Dx[n - j + 1, n - i + 1]
                end
            end
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
    Dxi = one(TR) ./ Dx  # elementwise reciprocal
    # Construct initial derivative matrix:
    D = Dw .* Dxi

    # Negative sum trick on the diagonal: D(i,i) = -sum of row i
    # This ensures each row sums to zero.
    for i in 1:n
        D[i, i] = -sum(D[i, :]) + D[i, i]  # add D[i,i] then subtract entire row
    end

    # Enforce a symmetry fix for even n on the diagonal entries:
    # (Matches the "Forcing symmetry for even n" from the original code.)
    # The original code effectively flips the sign for half of the diagonal indices for even n.
    # For simplicity, replicate that behavior exactly:
    halfN = div(n, 2)
    # Indices from end down to n-floor(n/2)+1:
    # That is from n down to n-halfN+1
    for idx in n:-1:(n - halfN + 1)
        D[idx, idx] = -D[idx, idx]
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
            D2[i, j] = 2.0 * D[i, j] * (dii - Dxi[i, j])
        end
    end
    # negative sum trick
    for i in 1:n
        D2[i, i] = -sum(D2[i, :]) + D2[i, i]
    end

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
        for i in 1:n
            D_next[i, i] = -sum(D_next[i, :]) + D_next[i, i]
        end
        D_current = D_next
    end

    return D_current
end
