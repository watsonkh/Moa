using StaticArrays

mutable struct Node
    position::MVector{3, Float64}
    displacement::MVector{3, Float64}
    velocity::MVector{3, Float64}
    @atomic force::MVector{3, Float64}
    volume::Float64
    material::Material
end

# Some constructors
Node(x::Float64, y::Float64, z::Float64, m::Material) = Node([x,y,z], [0,0,0], [0,0,0], [0,0,0], grid_spacing^3, m)
Node(pos::MVector{3, Float64}, m::Main.PD.Material, grid_spacing::Float64) = Node(pos, [0,0,0], [0,0,0], [0,0,0], grid_spacing^3, m)
