import Pkg
Pkg.activate(".")

using StaticArrays
import PyCall
import NearestNeighbors
using LinearAlgebra


# Global (for now) varibles
grid_spacing = 2.
horizon = grid_spacing * 3.014

mutable struct Material
    id::Int64
    density::Float64
    critical_strain::Float64
    bond_constant::Float64
end

mutable struct Node
    position::MVector{3, Float64}
    displacement::MVector{3, Float64}
    velocity::MVector{3, Float64}
    force::MVector{3, Float64}
    volume::Float64
    material::Material
end
# Some constructors
Node(x::Float64, y::Float64, z::Float64, m::Material) = Node([x,y,z], [0,0,0], [0,0,0], [0,0,0], grid_spacing^3, m)
Node(pos::MVector{3, Float64}, m::Material) = Node(pos, [0,0,0], [0,0,0], [0,0,0], grid_spacing^3, m)


mutable struct Bond
    from::Node
    to::Node
    isBroken::Bool
end

function apply_force(bond::Bond)
    initial_bond_length = norm(bond.to.position - bond.from.position)
    deformed_bond_vector = (bond.to.position + bond.to.displacement) - (bond.from.position + bond.from.displacement)
    deformed_bond_length = norm(deformed_bond_vector)
    strain = (deformed_bond_length - initial_bond_length) / initial_bond_length
    direction = deformed_bond_vector / deformed_bond_length

    bond.from.force = minimum([bond.from.material.bond_constant, bond.to.material.bond_constant]) * strain * direction * bond.to.volume * bond.from.volume
end

mutable struct Grid
    nodes::Vector{Node}
    bonds::Vector{Bond}
    material_default::Material
    materials::Vector{Material}
end




# Create Materials
materials = Vector{Material}()
push!(materials, Material(1, 1., 0.2, 75))
push!(materials, Material(2, 1., 0.2, 76))
println("Created ", length(materials), " materials!")

# Create Nodes (this process allows the nodes and position vector of vectors to point to the same thing)
positions = Vector{MVector{3, Float64}}()
push!(positions, [0.0, 0.0, 0.0])
push!(positions, [1.0, 0.0, 0.0])

nodes = Vector{Node}()
push!(nodes, Node(positions[1],materials[1]))
push!(nodes, Node(positions[2],materials[2]))
println("Created ", length(nodes), " nodes!")


# Create Bonds (double bonded)
bonds = Vector{Bond}()
tree = NearestNeighbors.BallTree(Vector{SVector{3, Float64}}(positions))
results = NearestNeighbors.inrange(tree, positions, horizon, false)
for (i, result) in enumerate(results)
    for neighbor in result
        if i != neighbor
            push!(bonds, Bond(nodes[i], nodes[neighbor], false))
        end
    end
end
println("Created ", length(bonds), " bonds!")

# Time iteration
node_dispalcement = Vector{Float64}() # used to record a node's dispalcement
nodes[2].displacement[1] = 0.05
nodes[1].displacement[1] = -0.05
# nodes[2].force[1] = 0.1
dt = grid_spacing / sqrt(maximum([material.bond_constant/material.density for material in materials])) * 0.01
for time_step in 1:1000000
    for node in nodes
        # Average velocity for the time step assuming constant acceleration
        node.velocity = node.velocity + (dt*0.5)*(node.force/(node.volume*node.material.density))

        # Update displacement with average velocity
        node.displacement = node.displacement + node.velocity * dt
    end
    # Apply displacement BCs
    # nodes[1].displacement[1] = 0.
    # Calculate forces
    for bond in bonds
        # Eventually replace node.force with atomic floats and use threading for this loop
        apply_force(bond)
    end
    for node in nodes
        # Calculate final velocity
        node.velocity = node.velocity + dt*0.5 * (node.force / (node.volume * node.material.density))
    end
    push!(node_dispalcement, nodes[2].displacement[1])
end

print(10000*dt)
using Plots
plot(node_dispalcement)
