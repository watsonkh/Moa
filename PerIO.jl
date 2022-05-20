
using TOML
import CSV
import NearestNeighbors

# include("Perijules.jl")
include("Bonds.jl")

input = TOML.parsefile("example_input.toml")



# nodes::Vector{Node} = []

# Read in all grids
for input_grid_object in input["Grid"]
    input_grid = CSV.File(input_grid_object["path"], stripwhitespace=true,  comment="#")

    # Required columns in grid file
    @assert :x in input_grid.names
    @assert :y in input_grid.names
    @assert :z in input_grid.names
    
    positions = convert(Matrix{Float64}, vcat(transpose(input_grid[:x]), transpose(input_grid[:y]), transpose(input_grid[:z])))
    tree = NearestNeighbors.BallTree(positions)
    result = NearestNeighbors.inrange(tree, positions, 3.014, false)
    println(typeof(result))
    n_points::Int64 = input_grid.rows
end

println("Ran successfully :)")