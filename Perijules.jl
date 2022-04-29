import Pkg
Pkg.activate(".")

import StaticArrays
import PyCall

mutable struct Material
    density::Float64
    critical_strain::Float64
    elastic_modulus::Float64
    grid_spacing::Float64
    horizon::Float64
end

mutable struct Node
    position::StaticArrays.SVector{3, Float64}
    displacement::StaticArrays.SVector{3, Float64}
    velocity::StaticArrays.SVector{3, Float64}
    force::StaticArrays.SVector{3, Float64}
    material::Material
end

mutable struct Bond
    from::Node
    to::Node
    isBroken::Bool
end

mutable struct Grid
    nodes::Vector{Node}
    bonds::Vector{Bond}
    material_default::Material
    materials::Vector{Material}
end

function parse_grid(filepath::String, skiprow::Bool)
    np = PyCall.pyimport("numpy")
    if skiprow
        input_matrix = np.loadtxt(filepath, skiprows=1)
    else
        input_matrix = np.loadtxt(filepath)
    end
    
    print(size(input_matrix))
end


function legacy_parse_grid(filepath::String)
    np = PyCall.pyimport("numpy")
    # Read in grid using numpy
    input_matrix = np.loadtxt(filepath, skiprows=1)

    # Initial positions are the first 3 columns in the matrix
    positions = input_matrix[:,1:3]
    displacements = zeros(size(positions))
    velocities = zeros(size(positions))
    forces = zeros(size(positions))

    # Material association numbers are last (4th) column in the matrix
    material_numbers = convert(Vector{Int32}, input_matrix[:,4])
    mat::Material = Material(1.,0.1,10000000,2.,2*3.015)
    nodes = []
    for i in 1:size(positions)[1]
        # n::Node = 
        append!(nodes,[Node(positions[i,:], displacements[i,:], velocities[i,:], forces[i,:], mat)])
    end

    return nodes
end




println("Data structure sizes:")
println("\tBond Size: ",        sizeof(Bond))
println("\tMaterial Size: ",    sizeof(Material))
println("\tNode Size: ",        sizeof(Node))
println("\tGrid Size: ",        sizeof(Grid), "\n\n")

# println("Nodes:\n",legacy_parse_grid("testgrid.grid"))
parse_grid("testgrid.grid", true)
# np = PyCall.pyimport("numpy")
# input_matrix = np.loadtxt("testgrid.grid", skiprows=1)

# println(typeof(arr))


# Read input file - materials, boundary conditions
    # Material
    # BC
        # Needs: Grid, Name(?), Type(?), capture criteria (inside volume, outside volume, crossing surface, etc.)
        # Types
            # Bond BC
                # Force Plane
            # Node BC
                # Displacement, Velocity, Acceleration, Staged Loading, etc
        # Pass it the grid and it handles everything?

# Read grid file - per-node properties (volume, position, displacement, nodal material, etc)




# Create materials

# Create nodes

# Create bonds

# Create boundary conditions


# Time iteration
