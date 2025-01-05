const AbstractFloatOrComplex = Union{AbstractFloat,Complex{<:AbstractFloat}}

function typeisfloat(::Type{T}) where {T<:Number}
    return T <: AbstractFloat
end

export typeisfloat, AbstractFloatOrComplex
