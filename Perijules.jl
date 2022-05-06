import Pkg
Pkg.activate(".")

import NearestNeighbors
using StaticArrays
using LinearAlgebra
import PyCall

# Global (for now) varibles
grid_spacing = 1.
horizon = grid_spacing * 3.014


module PD
    include("Materials.jl")
    include("Nodes.jl")
    include("Bonds.jl")
end


# TOOD for tensile bundle:
    # Nonlinear force response
    # Material interface
    # Staged loading BC
    # Contact algorithm




# Initialize grid

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
spatial = PyCall.pyimport_conda("scipy.spatial", "scipy", "anaconda")
tree = spatial.KDTree(positions);
py_neighbor_list = tree.query_ball_point(positions, horizon, workers=-1);
neighbor_list = Vector{Vector{Int64}}(py_neighbor_list);

for (i, result) in enumerate(neighbor_list)
    for neighbor in result
        n = neighbor+1
        if i != n
            push!(bonds, PD.Bond(nodes[i], nodes[n], false))
        end
    end
end
println("Created ", length(bonds), " bonds!")

## Plot initial grid
pv = PyCall.pyimport_conda("pyvista", "pyvista", "conda-forge")
point_cloud = pv.PolyData([node.position for node in nodes])
plotter = pv.Plotter()
plotter.add_mesh(point_cloud, scalars=material_numbers, cmap="plasma")
plotter.show()
##

SL_negative_points = [node for node in nodes if node.position[1] < -950]
SL_positive_points = [node for node in nodes if node.position[1] > 950]

FP_positive_bonds = [bond for bond in bonds if (bond.from.position[1] < 0 && bond.to.position[1] > 0)]
FP_negative_bonds = [bond for bond in bonds if (bond.from.position[1] > 0 && bond.to.position[1] < 0)]
FP_his = Vector{Float64}()

##

# Time iteration
ke_his = Vector{Float64}() # used to record kinetic energy

dt = grid_spacing / sqrt(maximum([material.bond_constant/material.density for material in materials])) * 0.01
for time_step in 1:1000
    Threads.@threads for node in nodes
        # Average velocity for the time step assuming constant acceleration
        node.velocity = node.velocity + (dt*0.5)*(node.force/(node.volume*node.material.density))

        # Update displacement with average velocity
        node.displacement += node.velocity * dt

        # Zero out force
        @atomic node.force = zeros(3)
    end

    # Calculate forces and break bonds
    Threads.@threads for bond in bonds
        if PD.should_break(bond)
            PD.break!(bond)
        end
        # Eventually replace node.force with atomic floats and use threading for this loop
        PD.apply_force(bond)
    end

    push!(FP_his, sum([PD.get_force(bond)[1] for bond in FP_positive_bonds]) -
     sum([PD.get_force(bond)[1] for bond in FP_negative_bonds]))

    Threads.@threads for node in nodes
        # Calculate final velocity
        node.velocity = node.velocity + dt*0.5 * (node.force / (node.volume * node.material.density))
    end

    kinetic_energy = sum([0.5 * node.volume * node.material.density * (norm(node.velocity)^2) for node in nodes])
    push!(ke_his, sum([0.5 * node.volume * node.material.density * (norm(node.velocity)^2) for node in nodes]))
    println(time_step)
end
println("Finished time loop!")

# Plot the two points' deformed position
import Plots
Plots.plot(FP_his, label=["Force plane apply_force"])