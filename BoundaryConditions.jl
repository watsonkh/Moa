module BoundaryConditions
using ..Materials
using ..Nodes
using ..Bonds

abstract type AbstractBoundaryCondition end

struct DisplacementBoundaryCondition <: AbstractBoundaryCondition
    nodes::Vector{Nodes.Node}
    displacement::Vector{Float64}
end

struct VelocityBoundaryCondition <: AbstractBoundaryCondition
    nodes::Vector{Nodes.Node}
    velocity::Vector{Float64}
end

# To add
#  Body force
#  (Surface pressure)

end