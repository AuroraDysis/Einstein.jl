const AbstractFloatOrComplex = Union{AbstractFloat,Complex{<:AbstractFloat}}

function typeisfloat(::Type{T}) where {T<:Number}
    return T <: AbstractFloat
end

function typetol(TR::Type{<:AbstractFloat})
    return real(oneunit(TR)) * (eps(real(one(TR))))^(4//5)
end

export typeisfloat, AbstractFloatOrComplex
export typetol
