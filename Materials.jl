module Materials

abstract type AbstractMaterial end

struct LinearEalstic <: AbstractMaterial
    id::Int64
    density::Float64
    critical_strain::Float64
    bond_constant::Float64
end

struct CustomMaterial <: AbstractMaterial
    id::Int64
    density::Float64
end


end
