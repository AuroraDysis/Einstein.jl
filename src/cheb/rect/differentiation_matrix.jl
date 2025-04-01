"""
    cheb_rect_differentiation_matrix_kind1([TF=Float64], m::Integer, n::Integer) where {TF<:AbstractFloat}
    cheb_rect_differentiation_matrix_kind1([TF=Float64], m::Integer, n::Integer, lower_bound::TF, upper_bound::TF) where {TF<:AbstractFloat}

Constructing a 1st-order rectangular differentiation matrix mapping from a 1st-kind grid to another 1st-kind grid.

# Arguments:
- `TF` : Type of the output matrix (default: Float64)
- `m` : Size of the output grid (number of rows).
- `n` : Size of the input grid (number of columns).

# References
- [chebfun/diffmat.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/diffmat.m)
"""
function cheb_rect_differentiation_matrix_kind1(::Type{TF}, m::Integer, n::Integer) where {TF<:AbstractFloat}
    # mapping-from grid (angles):
    T = angles(n)'        # Row vector of length n
    # difference between dimensions:
    c = n - m
    # mapping-to grid (angles):
    TAU = angles(TF, m)       # Column vector of length m

    # Denominator term
    denom = @. 2 * sin((T + TAU) / 2) * sin((TAU - T) / 2)

    # Sign matrix
    sgn = ones(TF, m, n)
    if isodd(c)
        sgn[1:2:m, 2:2:n] .= -1
        sgn[2:2:m, 1:2:n] .= -1
    else
        sgn[1:2:m, 1:2:n] .= -1
        sgn[2:2:m, 2:2:n] .= -1
    end

    # Matrix construction (fixed broadcasting)
    D = @. sgn * ((cos(c * TAU) / sin(TAU)) * sin(T) - sin(c * TAU) / n * sin(T) / denom) /
        denom

    # Negative-sum trick
    idx = [argmin(abs.(@view(denom[i, :]))) for i in 1:m]
    for i in 1:m
        D[i, idx[i]] = 0
        D[i, idx[i]] = -sum(@view D[i, :])
    end

    # Flipping trick
    if c == 1
        rot_D_180 = rot180(D)
        mask = rot180(tril(ones(Bool, m, n)))
        D[mask] .= -rot_D_180[mask]
    end

    return D
end

function cheb_rect_differentiation_matrix_kind1(m::Integer, n::Integer)
    return cheb_rect_differentiation_matrix_kind1(Float64, m, n)
end

function cheb_rect_differentiation_matrix_kind1(
    ::Type{TF}, m::Integer, n::Integer, lower_bound::TF, upper_bound::TF
) where {TF<:AbstractFloat}
    D = cheb_rect_differentiation_matrix_kind1(TF, m, n)
    D .*= 2 / (upper_bound - lower_bound)
    return D
end

function cheb_rect_differentiation_matrix_kind1(m::Integer, n::Integer, lower_bound::Float64, upper_bound::Float64)
    return cheb_rect_differentiation_matrix_kind1(Float64, m, n, lower_bound, upper_bound)
end

"""
    cheb_rect_differentiation_matrix_kind2([TF=Float64], m::Integer, n::Integer) where {TF<:AbstractFloat}
    cheb_rect_differentiation_matrix_kind2([TF=Float64], m::Integer, n::Integer, lower_bound::TF, upper_bound::TF) where {TF<:AbstractFloat}

Construct a 1st-order rectangular differentiation matrix mapping from a 2nd-kind grid to another 1st-kind grid.

# Arguments:
- `TF` : Type of the output matrix elements (e.g., Float64)
- `m` : Size of the output grid (number of rows)
- `n` : Size of the input grid (number of columns)

# References
- [chebfun/diffmat.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/diffmat.m)
"""
function cheb_rect_differentiation_matrix_kind2(::Type{TF}, m::Integer, n::Integer) where {TF<:AbstractFloat}
    nm1 = n - 1                     # For convenience
    cm1 = nm1 - m                   # Difference between dimensions
    t = points(TF, n)            # Second-kind grid
    tau = points(TF, m)          # First-kind grid
    T = angles(TF, n)         # Second-kind grid (angles)
    TAU = angles(TF, m)       # First-kind grid (angles)

    # Denominator term (explicit expression)
    denom = [2 * sin((t + tau) / 2) * sin((tau - t) / 2) for tau in TAU, t in T]

    # Numerator term
    numer = (1 .- tau * t') .* (cos.(cm1 .* TAU) ./ sin.(TAU))

    sgn = (-1)^cm1

    # Construct initial matrix
    if cm1 == 0
        D = numer ./ (denom .^ 2) ./ nm1
    else
        D = (sin.(cm1 .* TAU) .* ones(TF, 1, n)) ./ denom .+ numer ./ (denom .^ 2) ./ nm1
        D .*= sgn
    end

    # Scale first and last columns
    half = one(TF) / 2
    D[:, [1, n]] .*= half

    # Flipping trick for cm1 == 0
    if cm1 == 0
        ii = rot180(tril(ones(Bool, m, n)))
        rot90D = rot180(D)
        D[ii] .= sgn .* rot90D[ii]
    end

    # Sign adjustments
    D[1:2:end, 1:2:end] .*= -1
    D[2:2:end, 2:2:end] .*= -1

    # Negative sum trick
    idx = [argmin(abs.(denom[i, :])) for i in 1:m]
    for i in 1:m
        D[i, idx[i]] = 0
        D[i, idx[i]] = -sum(D[i, :])
    end

    # Fix corner values for cm1 == 0
    if cm1 == 0
        neg_quarter = -one(TF) / 4
        cornerVal = neg_quarter / (nm1 * sin(π / (2 * m)) * sin(π / (4 * m))^2)
        D[1, 1] = cornerVal
        D[end, end] = -cornerVal

        # Negative sum trick for corner entries
        if n > 2
            D[1, 2] = -sum(D[1, [1; 3:n]])
            D[end, end - 1] = -D[1, 2]
        end
    end

    return D
end

function cheb_rect_differentiation_matrix_kind2(m::Integer, n::Integer)
    return cheb_rect_differentiation_matrix_kind2(Float64, m, n)
end

function cheb_rect_differentiation_matrix_kind2(
    ::Type{TF}, m::Integer, n::Integer, lower_bound::TF, upper_bound::TF
) where {TF<:AbstractFloat}
    D = cheb_rect_differentiation_matrix_kind2(TF, m, n)
    D .*= 2 / (upper_bound - lower_bound)
    return D
end

function cheb_rect_differentiation_matrix_kind2(m::Integer, n::Integer, lower_bound::Float64, upper_bound::Float64)
    return cheb_rect_differentiation_matrix_kind2(Float64, m, n, lower_bound, upper_bound)
end

"""
    cheb_rect_differentiation_matrix([TF=Float64], m::Integer, n::Integer, p::Integer, kind::Integer) where {TF<:AbstractFloat}

Construct a p-th order rectangular differentiation matrix mapping from kind 1 or 2 Chebyshev grids to another 1st-kind grid.

# Arguments
- `TF` : Type of the output matrix elements (e.g., Float64)
- `m` : Size of the output grid (number of rows)
- `n` : Size of the input grid (number of columns)
- `p` : Order of differentiation
- `kind` : Kind of Chebyshev grid (1 or 2)
"""
function cheb_rect_differentiation_matrix(
    ::Type{TF}, m::Integer, n::Integer, p::Integer, kind::Integer
)::Matrix{TF} where {TF<:AbstractFloat}
    @argcheck p >= 2 "p must be at least 2"
    @argcheck kind == 1 || kind == 2 "kind must be 1 or 2"

    # Initialize sign vector
    sgn = ones(TF, n)
    sgn[1:2:n] .= -1

    if kind == 1
        # First-kind grid
        T = angles(TF, n)
        D = cheb_rect_differentiation_matrix_kind1(TF, m, n)
        a = vcat(zeros(TF, n), one(TF))
        sgn_coeff = (-1)^(n - 1) / TF(n)
        @. sgn = sgn_coeff * sgn * sin(T)
    else
        # Second-kind grid
        T = angles(TF, n)
        D = cheb_rect_differentiation_matrix_kind2(TF, m, n)
        a = vcat(zeros(TF, n - 2), -one(TF), zero(TF), one(TF))
        sgn .*= (-1)^(n - 1) / TF(2 * (n - 1))
        sgn[[1, n]] .*= one(TF) / 2
    end

    # Setup grids
    tau = points(TF, m)
    TAU = angles(TF, m)
    a = cheb_coeffs_diff(a)

    # Compute denominator matrix
    denom = [2 * sin((tau + t) / 2) * sin((tau - t) / 2) for tau in TAU, t in T]
    idx = @inbounds [argmin(abs.(@view(denom[i, :]))) for i in 1:m]

    # Higher-order derivatives
    for l in 2:p
        a = cheb_coeffs_diff(a)
        Tt = cheb_coeffs_eval(a, tau)
        D .= (Tt .* sgn' .+ l .* D) ./ denom
    end

    # Apply negative-sum trick
    @inbounds for i in 1:m
        row_sum = sum(@view(D[i, :]))
        D[i, idx[i]] = -row_sum + D[i, idx[i]]
    end

    return D
end

function cheb_rect_differentiation_matrix(m::Integer, n::Integer, p::Integer, kind::Integer)
    return cheb_rect_differentiation_matrix(Float64, m, n, p, kind)
end

function cheb_rect_differentiation_matrix(
    ::Type{TF}, m::Integer, n::Integer, p::Integer, kind::Integer, lower_bound::TF, upper_bound::TF
)::Matrix{TF} where {TF<:AbstractFloat}
    D = cheb_rect_differentiation_matrix(TF, m, n, p, kind)
    scale = (2 / (upper_bound - lower_bound))^p
    D .*= scale
    return D
end

function cheb_rect_differentiation_matrix(
    m::Integer, n::Integer, p::Integer, kind::Integer, lower_bound::Float64, upper_bound::Float64
)
    return cheb_rect_differentiation_matrix(Float64, m, n, p, kind, lower_bound, upper_bound)
end

export cheb_rect_differentiation_matrix_kind1, cheb_rect_differentiation_matrix_kind2, cheb_rect_differentiation_matrix
