"""
    cheb_rectdiff1([TR=Float64], m::Integer, n::Integer) where {TR<:AbstractFloat}
    cheb_rectdiff1([TR=Float64], m::Integer, n::Integer, x_min::TR, x_max::TR) where {TR<:AbstractFloat}

Constructing a 1st-order rectangular differentiation matrix mapping from a 1st-kind grid

# Arguments:
- `m` : Size of the output grid (number of rows).
- `n` : Size of the input grid (number of columns).

# References
- [chebfun/diffmat.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/diffmat.m)
"""
function cheb_rectdiff1(::Type{TR}, m::Integer, n::Integer) where {TR<:AbstractFloat}
    # mapping-from grid (angles):
    T = cheb1_angles(n)'        # Row vector of length n
    # difference between dimensions:
    c = n - m
    # mapping-to grid (angles):
    TAU = cheb1_angles(TR, m)       # Column vector of length m

    # Denominator term
    denom = @. 2 * sin((T + TAU) / 2) * sin((TAU - T) / 2)

    # Sign matrix
    sgn = ones(TR, m, n)
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

function cheb_rectdiff1(m::Integer, n::Integer)
    return cheb_rectdiff1(Float64, m, n)
end

function cheb_rectdiff1(
    ::Type{TR}, m::Integer, n::Integer, x_min::TR, x_max::TR
) where {TR<:AbstractFloat}
    D = cheb_rectdiff1(TR, m, n)
    D .*= 2 / (x_max - x_min)
    return D
end

function cheb_rectdiff1(m::Integer, n::Integer, x_min::Float64, x_max::Float64)
    return cheb_rectdiff1(Float64, m, n, x_min, x_max)
end

"""
    cheb_rectdiff2([TR=Float64], m::Integer, n::Integer) where {TR<:AbstractFloat}
    cheb_rectdiff2([TR=Float64], m::Integer, n::Integer, x_min::TR, x_max::TR) where {TR<:AbstractFloat}

Construct a 1st-order rectangular differentiation matrix mapping from a 2nd-kind grid.

# Arguments:
- `m` : Size of the output grid (number of rows)
- `n` : Size of the input grid (number of columns)

# References
- [chebfun/diffmat.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/diffmat.m)
"""
function cheb_rectdiff2(::Type{TR}, m::Integer, n::Integer) where {TR<:AbstractFloat}
    nm1 = n - 1                     # For convenience
    cm1 = nm1 - m                   # Difference between dimensions
    t = cheb2_points(TR, n)            # Second-kind grid
    tau = cheb1_points(TR, m)          # First-kind grid
    T = cheb2_angles(TR, n)         # Second-kind grid (angles)
    TAU = cheb1_angles(TR, m)       # First-kind grid (angles)

    # Denominator term (explicit expression)
    denom = [2 * sin((t + tau) / 2) * sin((tau - t) / 2) for tau in TAU, t in T]

    # Numerator term
    numer = (1 .- tau * t') .* (cos.(cm1 .* TAU) ./ sin.(TAU))

    sgn = (-1)^cm1

    # Construct initial matrix
    if cm1 == 0
        D = numer ./ (denom .^ 2) ./ nm1
    else
        D = (sin.(cm1 .* TAU) .* ones(TR, 1, n)) ./ denom .+ numer ./ (denom .^ 2) ./ nm1
        D .*= sgn
    end

    # Scale first and last columns
    half = one(TR) / 2
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
        neg_quarter = -one(TR) / 4
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

function cheb_rectdiff2(m::Integer, n::Integer)
    return cheb_rectdiff2(Float64, m, n)
end

function cheb_rectdiff2(
    ::Type{TR}, m::Integer, n::Integer, x_min::TR, x_max::TR
) where {TR<:AbstractFloat}
    D = cheb_rectdiff2(TR, m, n)
    D .*= 2 / (x_max - x_min)
    return D
end

function cheb_rectdiff2(m::Integer, n::Integer, x_min::Float64, x_max::Float64)
    return cheb_rectdiff2(Float64, m, n, x_min, x_max)
end

"""
    cheb_rectdiff([TR=Float64], m::Integer, n::Integer, p::Integer, kind::Integer) where {TR<:AbstractFloat}

Construct a p-th order rectangular differentiation matrix mapping between Chebyshev grids.

# Arguments
- `m` : Size of the output grid (number of rows)
- `n` : Size of the input grid (number of columns)
- `p` : Order of differentiation
- `kind` : Kind of Chebyshev grid (1 or 2)
"""
function cheb_rectdiff_rec(
    ::Type{TR}, m::Integer, n::Integer, p::Integer, kind::Integer
)::Matrix{TR} where {TR<:AbstractFloat}
    @argcheck p >= 2 "p must be at least 2"
    @argcheck kind == 1 || kind == 2 "kind must be 1 or 2"

    # Initialize sign vector
    sgn = ones(TR, n)
    sgn[1:2:n] .= -1

    if kind == 1
        # First-kind grid
        T = cheb1_angles(TR, n)
        D = cheb_rectdiff1(TR, m, n)
        a = vcat(zeros(TR, n), one(TR))
        sgn_coeff = (-1)^(n - 1) / TR(n)
        @. sgn = sgn_coeff * sgn * sin(T)
    else
        # Second-kind grid
        T = cheb2_angles(TR, n)
        D = cheb_rectdiff2(TR, m, n)
        a = vcat(zeros(TR, n - 2), -one(TR), zero(TR), one(TR))
        sgn .*= (-1)^(n - 1) / TR(2 * (n - 1))
        sgn[[1, n]] .*= one(TR) / 2
    end

    # Setup grids
    tau = cheb1_points(TR, m)
    TAU = cheb1_angles(TR, m)
    a = cheb_diff(a)

    # Compute denominator matrix
    denom = [2 * sin((tau + t) / 2) * sin((tau - t) / 2) for tau in TAU, t in T]
    idx = @inbounds [argmin(abs.(@view(denom[i, :]))) for i in 1:m]

    # Higher-order derivatives
    for l in 2:p
        a = cheb_diff(a)
        Tt = cheb_clenshaw(a, tau)
        D .= (Tt .* sgn' .+ l .* D) ./ denom
    end

    # Apply negative-sum trick
    @inbounds for i in 1:m
        row_sum = sum(@view(D[i, :]))
        D[i, idx[i]] = -row_sum + D[i, idx[i]]
    end

    return D
end

function cheb_rectdiff_rec(m::Integer, n::Integer, p::Integer, kind::Integer)
    return cheb_rectdiff_rec(Float64, m, n, p, kind)
end

function cheb_rectdiff_rec(
    ::Type{TR}, m::Integer, n::Integer, p::Integer, kind::Integer, x_min::TR, x_max::TR
)::Matrix{TR} where {TR<:AbstractFloat}
    D = cheb_rectdiff_rec(TR, m, n, p, kind)
    scale = (2 / (x_max - x_min))^p
    D .*= scale
    return D
end

function cheb_rectdiff_rec(
    m::Integer, n::Integer, p::Integer, kind::Integer, x_min::Float64, x_max::Float64
)
    return cheb_rectdiff_rec(Float64, m, n, p, kind, x_min, x_max)
end

export cheb_rectdiff1, cheb_rectdiff2, cheb_rectdiff_rec
