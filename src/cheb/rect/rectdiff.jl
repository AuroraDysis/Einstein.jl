"""
    cheb_rectdiff1([TR=Float64], m::TI, n::TI) where {TR<:AbstractFloat,TI<:Integer}

Constructing a 1st-order rectangular differentiation matrix mapping from a 1st-kind grid

Arguments:
- `m` : Size of the output grid (number of rows).
- `n` : Size of the input grid (number of columns).

Returns:
- `D` : The rectangular differentiation matrix (size `m x n`).
"""
function cheb_rectdiff1(::Type{TR}, m::TI, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    # mapping-from grid (angles):
    T = cheb1_angles(n)'        # Row vector of length n
    # difference between dimensions:
    c = n - m
    # mapping-to grid (angles):
    TAU = cheb1_angles(m)       # Column vector of length m

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

"""
    cheb_rectdiff2([TR=Float64], m::TI, n::TI) where {TR<:AbstractFloat,TI<:Integer}

Construct a 1st-order rectangular differentiation matrix mapping from a 2nd-kind grid.

Arguments:
- `m` : Size of the output grid (number of rows)
- `n` : Size of the input grid (number of columns)

Returns:
- `D` : The rectangular differentiation matrix (size `m x n`)
"""
function cheb_rectdiff2(::Type{TR}, m::TI, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    nm1 = n - 1                     # For convenience
    cm1 = nm1 - m                   # Difference between dimensions
    t = cheb2_pts(TR, n)            # Second-kind grid
    tau = cheb1_pts(TR, m)          # First-kind grid
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

function cheb_rectdiff1(m::TI, n::TI) where {TI<:Integer}
    return cheb_rectdiff1(Float64, m, n)
end

function cheb_rectdiff2(m::TI, n::TI) where {TI<:Integer}
    return cheb_rectdiff2(Float64, m, n)
end

export cheb_rectdiff1, cheb_rectdiff2

@testset "cheb_rectdiff1" begin
    tol = 100 * eps(Float64)

    D45 = cheb_rectdiff1(Float64, 4, 5)
    D45_res = [
        -4.57116486120251 6.35944256641146 -2.70279423558635 1.31914227171008 -0.404625741332680
        -0.0484519638509037 -1.91971151095453 2.34412101225198 -0.504026090801922 0.128068553355378
        -0.128068553355378 0.504026090801922 -2.34412101225198 1.91971151095453 0.0484519638509037
        0.404625741332680 -1.31914227171008 2.70279423558635 -6.35944256641146 4.57116486120251
    ]
    @test isapprox(D45, D45_res; atol=tol)

    D56 = cheb_rectdiff1(Float64, 5, 6)
    D56_res = [
        -6.71910606480245 9.53284806144464 -4.39831637342411 2.49110926664733 -1.32569196231062 0.419157072445230
        -0.253223738073512 -2.39105495582767 3.33681201184053 -1.01065674661971 0.453608214480045 -0.135484785799679
        -0.0462335694140099 0.235702260395516 -2.40325617336917 2.40325617336917 -0.235702260395516 0.0462335694140099
        0.135484785799679 -0.453608214480045 1.01065674661971 -3.33681201184053 2.39105495582767 0.253223738073512
        -0.419157072445230 1.32569196231062 -2.49110926664733 4.39831637342411 -9.53284806144464 6.71910606480245
    ]
    @test isapprox(D56, D56_res; atol=tol)
end

@testset "cheb_rectdiff2" begin
    tol = 1000 * eps(Float64)

    D45 = cheb_rectdiff2(Float64, 4, 5)
    D45_res = [
        -4.29110266916749 4.82023271093930 -0.765366864730180 0.406019148566206 -0.169782325607842
        0.219172839560929 -1.87528541910585 1.84775906502257 -0.289498981478942 0.0978524960012859
        -0.0978524960012859 0.289498981478942 -1.84775906502257 1.87528541910585 -0.219172839560929
        0.169782325607842 -0.406019148566206 0.765366864730180 -4.82023271093930 4.29110266916749
    ]
    @test isapprox(D45, D45_res; atol=tol)

    D56 = cheb_rectdiff2(Float64, 5, 6)
    D56_res = [
        -6.61184642477616 7.39689062697822 -1.10865099995204 0.527416977547751 -0.369672519021278 0.165862339223506
        0.299860202570142 -2.64910740362000 2.60334916729954 -0.363213978699957 0.186960573879793 -0.0778485614295215
        -0.100000000000000 0.305572809000084 -2.09442719099992 2.09442719099992 -0.305572809000084 0.100000000000000
        0.0778485614295215 -0.186960573879793 0.363213978699957 -2.60334916729954 2.64910740362000 -0.299860202570142
        -0.165862339223506 0.369672519021278 -0.527416977547751 1.10865099995204 -7.39689062697822 6.61184642477616
    ]
    @test isapprox(D56, D56_res; atol=tol)
end
