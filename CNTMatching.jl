module PD
    include("Materials.jl")
    include("Nodes.jl")
    include("CNTBonds.jl")
    include("CellList.jl")
end

import Pkg
Pkg.activate(".")

import NearestNeighbors
using StaticArrays
using LinearAlgebra
import PyCall


# Global (for now) varibles
grid_spacing = 2.
horizon = grid_spacing * 3.014

# Create material
materials = Vector{PD.Material}()
push!(materials, PD.Material(1, 1., 0.201, 0.48828125/75.0))

# Create nodes
nodes = Vector{PD.Node}()
for x in 0:2:10
    println(x)
    push!(nodes, PD.Node(Float64(x), 0., 0., materials[1], grid_spacing))
end

# Create Bonds (double bonded)
cntBonds = Vector{PD.CNTBond}()
# interfaceBonds = Vector{PD.CNTBond}()

cell_list = PD.create_cell_list(nodes, horizon)
for node in nodes
    for other in PD.sample_cell_list(cell_list, node, horizon)
        if node.material.id == other.material.id
            push!(cntBonds, PD.CNTBond(node, other, false))
        end
    end
end

# Force Plane
FP_positive_bonds = [bond for bond in cntBonds if (bond.from.position[1] < 5 && bond.to.position[1] > 5)]
FP_negative_bonds = [bond for bond in cntBonds if (bond.from.position[1] > 5 && bond.to.position[1] < 5)]
FP_his = Vector{Float64}()
disp_his = Vector{Float64}()

SL_increment = 0.01
dt = grid_spacing / sqrt(maximum([material.bond_constant/material.density for material in materials])) * 0.01

for time_step in 1:300

    Threads.@threads for node in nodes
        # Average velocity for the time step assuming constant acceleration
        node.velocity = node.velocity + (dt*0.5)*(node.force/(node.volume*node.material.density))

        # Update displacement with average velocity
        node.displacement += node.velocity * dt

        # Zero out force
        @atomic node.force = zeros(3)
    end

    for node in nodes
        node.displacement = node.position * time_step * SL_increment * 0.1
        # Zero velocity (known static equlibrium)
        node.velocity = zeros(3)
    end

    # Break bonds
    Threads.@threads for bond in cntBonds
        if PD.should_break(bond)
            PD.break!(bond)
        end
    end

    # Apply bond force to nodes
    Threads.@threads for bond in cntBonds
        PD.apply_force(bond)
    end

    # Calculate final velocity
    Threads.@threads for node in nodes
        node.velocity = node.velocity + dt*0.5 * (node.force / (node.volume * node.material.density))
    end

    push!(FP_his, (sum([PD.get_force(bond)[1] for bond in FP_positive_bonds]) - sum([PD.get_force(bond)[1] for bond in FP_negative_bonds]))*0.5)
    push!(disp_his, time_step * SL_increment)
end

import Plots
Plots.xlabel!("Strain [nm/nm]")
Plots.ylabel!("Force [nm]")
print("Max force = ", maximum(FP_his))

Plots.plot(disp_his/10.0, FP_his, label="Peridynamic force response")
Plots.plot!([row[1] for row in PD.strainForce], [row[2] for row in PD.strainForce],label="MD Simulation")