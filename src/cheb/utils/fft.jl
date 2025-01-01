function plan_fft_measure!(tmp::Vector{Complex{T}}) where {T}
    plan = if T <: FFTW.fftwReal
        plan_fft!(tmp; flags=FFTW.MEASURE)
    else
        plan_fft!(tmp)
    end
    return plan
end

function plan_ifft_measure!(tmp::Vector{Complex{T}}) where {T}
    plan = if T <: FFTW.fftwReal
        plan_ifft!(tmp; flags=FFTW.MEASURE)
    else
        plan_ifft!(tmp)
    end
    return plan
end
