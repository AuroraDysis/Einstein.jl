function typetol(TF::Type{<:AbstractFloat})
    return real(oneunit(TF)) * (eps(real(one(TF))))^(4//5)
end

export typetol
