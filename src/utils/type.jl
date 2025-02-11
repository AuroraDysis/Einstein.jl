const AbstractFloatOrComplex = Union{AbstractFloat,Complex{<:AbstractFloat}}

function typeisfloat(::Type{T}) where {T<:Number}
    return T <: AbstractFloat
end

function typetol(TF::Type{<:AbstractFloat})
    return real(oneunit(TF)) * (eps(real(one(TF))))^(4//5)
end

export typeisfloat, AbstractFloatOrComplex
export typetol
