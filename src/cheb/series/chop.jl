function _cheb_series_chop_tol_impl!(
    coeffs::AbstractVector{TFC}, tol::TF
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}}
    @inbounds for k in length(coeffs):-1:1
        if abs(coeffs[k]) > tol
            return k
        end
    end
    return 1
end

function _cheb_series_chop_impl!(
    coeffs::AbstractVector{TFC}, tol::TF, envelope::AbstractVector{TF}
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}}
    n = length(coeffs)

    # Step 1: Convert COEFFS to a new monotonically nonincreasing vector
    # ENVELOPE normalized to begin with the value 1.

    @.. envelope = abs(coeffs)

    # Compute the accumulated maximum of the absolute values of the coefficients
    @inbounds for j in (n - 1):-1:1
        envelope[j] = max(envelope[j], envelope[j + 1])
    end

    if iszero(envelope[1])
        return 1
    end

    @.. envelope /= envelope[1]

    # Step 2: Scan ENVELOPE for a value PLATEAUPOINT, the first point J-1, if any,
    # that is followed by a plateau.  A plateau is a stretch of coefficients
    # ENVELOPE(J),...,ENVELOPE(J2), J2 = round(1.25*J+5) <= N, with the property
    # that ENVELOPE(J2)/ENVELOPE(J) > R.  The number R ranges from R = 0 if
    # ENVELOPE(J) = TOL up to R = 1 if ENVELOPE(J) = TOL^(2/3).  Thus a potential
    # plateau whose starting value is ENVELOPE(J) ~ TOL^(2/3) has to be perfectly
    # flat to count, whereas with ENVELOPE(J) ~ TOL it doesn't have to be flat at
    # all.  If a plateau point is found, then we know we are going to chop the
    # vector, but the precise chopping point CUTOFF still remains to be determined
    # in Step 3.

    plateau_point = 0

    for j in 2:n
        j2 = round(Int, 1.25 * j + 5)

        if j2 > n
            # there is no plateau: exit
            return n
        end

        e1 = envelope[j]
        e2 = envelope[j2]

        r = 3 * (1 - log(e1) / log(tol))
        plateau = iszero(e1) || (e2 / e1 > r)

        if plateau
            # a plateau has been found: go to Step 3
            plateau_point = j - 1
            break
        end
    end

    # Step 3: fix CUTOFF at a point where ENVELOPE, plus a linear function
    # included to bias the result towards the left end, is minimal.
    #
    # Some explanation is needed here.  One might imagine that if a plateau is
    # found, then one should simply set CUTOFF = PLATEAUPOINT and be done, without
    # the need for a Step 3. However, sometimes CUTOFF should be smaller or larger
    # than PLATEAUPOINT, and that is what Step 3 achieves.
    #
    # CUTOFF should be smaller than PLATEAUPOINT if the last few coefficients made
    # negligible improvement but just managed to bring the vector ENVELOPE below the
    # level TOL^(2/3), above which no plateau will ever be detected.  This part of
    # the code is important for avoiding situations where a coefficient vector is
    # chopped at a point that looks "obviously wrong" with PLOTCOEFFS.
    #
    # CUTOFF should be larger than PLATEAUPOINT if, although a plateau has been
    # found, one can nevertheless reduce the amplitude of the coefficients a good
    # deal further by taking more of them.  This will happen most often when a
    # plateau is detected at an amplitude close to TOL, because in this case, the
    # "plateau" need not be very flat.  This part of the code is important to
    # getting an extra digit or two beyond the minimal prescribed accuracy when it
    # is easy to do so.

    if plateau_point != 0 && iszero(envelope[plateau_point])
        return plateau_point
    end

    tol_7_6 = tol^(7//6)

    j3 = sum(envelope .≥ tol_7_6)
    if j3 < j2
        j2 = j3 + 1
        envelope[j2] = tol_7_6
    end

    envelope_j2 = @view(envelope[1:j2])
    range_j2 = range(0; stop=-log10(tol) / 3, length=j2)

    @.. envelope_j2 = log10(envelope_j2) + range_j2
    d = argmin(envelope_j2)
    return max(d - 1, 1)
end

"""
    cheb_series_chop(coeffs::AbstractVector{TF}, tol::TF = eps(TF)) where {TF<:AbstractFloat}

Determine a suitable cutoff index for a coefficient vector using the "standard" chopping rule [aurentz2015choppingchebyshevseries](@cite).

# References
- [aurentz2015choppingchebyshevseries](@citet*) Aurentz, J. L., & Trefethen, L. N. (2015). Chopping a Chebyshev series. ACM Transactions on Mathematical Software (TOMS), 41(4), 1-18.
- [chebfun/standardChop.m](https://github.com/chebfun/chebfun/blob/master/%40chebtech/standardChop.m) (Implementation in Chebfun, MATLAB)
"""
function cheb_series_chop(
    coeffs::AbstractVector{TFC},
    tol::TF=eps(TF);
    envelope::AbstractVector{TF}=Vector{TF}(undef, length(coeffs)),
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}}
    @boundscheck begin
        @argcheck !isempty(coeffs) "coeffs must not be empty"
        @argcheck 0 < tol < 1 "tol must be between 0 and 1 (exclusive)"
        @argcheck length(envelope) >= length(coeffs) "envelope must be at least as long as coeffs"
    end

    n = length(coeffs)

    if n < 17
        # If the coeffs is too short, use the tolerance-based implementation.
        return _cheb_series_chop_tol_impl!(coeffs, tol)
    end

    envelope_view = if length(envelope) > n
        @view(envelope[1:n])
    else
        envelope
    end

    return _cheb_series_chop_impl!(coeffs, tol, envelope_view)
end

export cheb_series_chop
