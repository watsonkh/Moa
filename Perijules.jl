    
module PD
    include("Materials.jl")
    include("Nodes.jl")
    include("Bonds.jl")
    include("CellList.jl")
end

import Pkg
Pkg.activate(".")
Pkg.add("StaticArrays")


import NearestNeighbors
using StaticArrays
using LinearAlgebra
import PyCall


# Global (for now) varibles
grid_spacing = 2.
horizon = grid_spacing * 3.014


# TOOD for tensile bundle:
    # Nonlinear force response
    # Material interface
    # Contact algorithm


# Read in input grid
np = PyCall.pyimport("numpy")
data = np.loadtxt("testgrid.grid", skiprows=1)

# Create materials
material_numbers = Vector{Int64}(data[:,4])
unique_material_numbers = sort(unique!(Vector{Int64}(data[:,4])))
materials = Vector{PD.Material}()
for material_number in unique_material_numbers
    push!(materials, PD.Material(material_number, 1., 0.1, 75))
end

# Create Nodes (this process allows the nodes and position vector of vectors to point to the same thing)
positions = data[:,1:3]
nodes = Vector{PD.Node}()
for (i, position) in enumerate(eachrow(positions))
    material = [material for material in materials if material.id == material_numbers[i]][1]
    push!(nodes, PD.Node(MVector{3, Float64}(position), material, grid_spacing))
end
println("Created ", length(nodes), " nodes!")

# Create Bonds (double bonded)
bonds = Vector{PD.Bond}()
@time cell_list = PD.create_cell_list(nodes, horizon);
@time cell_list = PD.create_cell_list_new(nodes, horizon);
for node in nodes
    # Can have lots of optimizations here
    for other in PD.sample_cell_list(cell_list, node, horizon)
        push!(bonds, PD.Bond(node, other, false))
    end
end
println("Created ", length(bonds), " bonds!")


# Staged Loading

# 1 - displace and relax
# 2 - fail and relax
SL_stage = 1
SL_positive_points = [node for node in nodes if node.position[1] >  950]
SL_negative_points = [node for node in nodes if node.position[1] < -950]
SL_iteration_count = 1
SL_increment = 10.
SL_last_dynamic_time_step = 1

# Force Plane
FP_positive_bonds = [bond for bond in bonds if (bond.from.position[1] < 0 && bond.to.position[1] > 0)]
FP_negative_bonds = [bond for bond in bonds if (bond.from.position[1] > 0 && bond.to.position[1] < 0)]
FP_his = Vector{Float64}()
disp_his = Vector{Float64}()
# Time iteration
ke_his = Vector{Float64}() # used to record kinetic energy

dt = grid_spacing / sqrt(maximum([material.bond_constant/material.density for material in materials])) * 0.01


for time_step in 1:8000

    Threads.@threads for node in nodes
        # Average velocity for the time step assuming constant acceleration
        node.velocity = node.velocity + (dt*0.5)*(node.force/(node.volume*node.material.density))

        # Update displacement with average velocity
        node.displacement += node.velocity * dt

        # Zero out force
        @atomic node.force = zeros(3)
    end

    # Staged loading: constrain tabs
    Threads.@threads for node in SL_positive_points
        node.displacement = [1.,0.,0.] * SL_iteration_count * SL_increment
        node.velocity = zeros(3)
    end
    Threads.@threads for node in SL_negative_points
        node.displacement = -[1.,0.,0.] * SL_iteration_count * SL_increment
        node.velocity = zeros(3)
    end

    if SL_stage == 2
        # Break bonds
        Threads.@threads for bond in bonds
            if PD.should_break(bond)
                PD.break!(bond)
            end
        end
    end

    # Apply bond force to nodes
    Threads.@threads for bond in bonds
        PD.apply_force(bond)
    end


    Threads.@threads for node in nodes
        # Calculate final velocity
        node.velocity = node.velocity + dt*0.5 * (node.force / (node.volume * node.material.density))
        node.velocity *= 0.99
    end

    kinetic_energy = sum([0.5 * node.volume * node.material.density * (norm(node.velocity)^2) for node in nodes])
    if kinetic_energy < 0.0025
        if SL_stage == 1 time_step - SL_last_dynamic_time_step > 50
            SL_stage = 2
            SL_last_dynamic_time_step = time_step
        end
        if SL_stage == 2 && time_step - SL_last_dynamic_time_step > 50
            SL_stage = 1
            SL_last_dynamic_time_step = time_step
            SL_iteration_count += 1
            println("NEW STAGE ##############################")
            push!(FP_his, (sum([PD.get_force(bond)[1] for bond in FP_positive_bonds]) - sum([PD.get_force(bond)[1] for bond in FP_negative_bonds]))*0.5)
            push!(disp_his, SL_iteration_count * SL_increment)
        end
    end
    push!(ke_his, kinetic_energy)

    println(time_step, ": ", kinetic_energy)
end
# println("Finished time loop!")

# Plot the two points' deformed position
import Plots
# Plots.plot(FP_his, label=["Force plane apply_force"])
Plots.plot(ke_his, label="Kinetic Energy")

# SL_iteration_count


## Plotting with pyvista
pv = PyCall.pyimport_conda("pyvista", "pyvista", "conda-forge")
point_cloud = pv.PolyData([node.position+node.displacement for node in nodes])
plotter = pv.Plotter()
plotter.add_mesh(point_cloud, scalars=[node.displacement[1] for node in nodes], cmap="plasma")
plotter.show()
