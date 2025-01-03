function plan_fft_measure!(tmp::Vector{Complex{T}}) where {T<:FFTW.fftwReal}
    return plan_fft!(tmp; flags=FFTW.MEASURE)
end

function plan_fft_measure!(tmp::Vector{Complex{T}}) where {T<:AbstractFloat}
    return plan_fft!(tmp)
end

function plan_ifft_measure!(tmp::Vector{Complex{T}}) where {T<:FFTW.fftwReal}
    return plan_ifft!(tmp; flags=FFTW.MEASURE)
end

function plan_ifft_measure!(tmp::Vector{Complex{T}}) where {T<:AbstractFloat}
    return plan_ifft!(tmp)
end
