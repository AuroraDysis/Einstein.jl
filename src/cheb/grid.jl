@enum ChebGridKind begin
    FirstKind = 1
    SecondKind = 2
end

function cheb_grid_first_kind(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    grid = zeros(TR, n)

    pi_over_2n = convert(TR, pi) / (2 * n)
    @inbounds begin
        for k in 0:(n - 1)
            grid[k + 1] = -cos((2 * k + 1) * pi_over_2n)
        end
    end

    return grid
end

function cheb_grid_second_kind(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    grid = zeros(TR, n)

    pi_over_nm1 = convert(TR, pi) / (n - 1)
    @inbounds begin
        grid[1] = -1
        grid[n] = 1
        for k in 1:(n - 2)
            grid[k + 1] = -cos(k * pi_over_nm1)
        end
    end

    return grid
end

function cheb_grid(
    ::Type{TR}, n::TI, kind::ChebGridKind=SecondKind
) where {TR<:AbstractFloat,TI<:Integer}
    if kind == FirstKind
        return cheb_grid_first_kind(TR, n)
    elseif kind == SecondKind
        return cheb_grid_second_kind(TR, n)
    else
        throw(ArgumentError("kind must be FirstKind or SecondKind"))
    end
end

export cheb_grid
export ChebGridKind, FirstKind, SecondKind
