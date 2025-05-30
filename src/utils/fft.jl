function plan_fft_measure!(tmp::Vector{Complex{TF}}) where {TF<:AbstractFloat}
    @argcheck length(tmp) > 1 "tmp must have at least two elements"

    fill!(tmp, zero(Complex{TF}))

    if TF <: FFTW.fftwReal
        return plan_fft!(tmp; flags=FFTW.MEASURE)
    else
        return plan_fft!(tmp)
    end
end

function plan_ifft_measure!(tmp::Vector{Complex{TF}}) where {TF<:AbstractFloat}
    @argcheck length(tmp) > 1 "tmp must have at least two elements"

    fill!(tmp, zero(Complex{TF}))

    if TF <: FFTW.fftwReal
        return plan_ifft!(tmp; flags=FFTW.MEASURE)
    else
        return plan_ifft!(tmp)
    end
end

export plan_fft_measure!, plan_ifft_measure!
