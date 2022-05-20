using TOML
using CSV

function parse_material(inputDict)
    # A material needs at least these 3 properties
    @assert haskey(inputDict, "type")
    @assert haskey(inputDict, "density")
    @assert haskey(inputDict, "id")

    if inputDict["type"] == "LinearElastic"
        return Materials.LinearEalstic(inputDict["id"], inputDict["density"], inputDict["critical_strain"], inputDict["bond_constant"])
    elseif inputDict["type"] == "CustomMaterial"
        return Materials.CustomMaterial(inputDict["id"], inputDict["density"])
    else
        # Material type not known
        println("#### Unknown type from material: ", inputDict["id"])
        throw(Exception)
    end
end


function parse_input(path::String)
    input = TOML.parsefile(path)
    
    println("Parsing input...")

    # Parse global scalars
    @show global gridspacing = input["gridspacing"]
    @show global horizon = input["horizon"]

    # Parse materials
    materials::Vector{Materials.AbstractMaterial} = Vector{Materials.AbstractMaterial}()

    defaultMaterial = parse_material(input["MaterialDefault"])
    push!(materials, defaultMaterial)

    for materialInput in input["Material"]
        mat::Materials.AbstractMaterial = parse_material(materialInput)
        # Each material should have a unique id
        @assert !(mat.id in [material.id for material in materials])
        push!(materials, mat)
    end

    # Parse grids
    nodes::Vector{Nodes.AbstractNode} = Vector{Nodes.AbstractNode}()
    for input_grid_object in input["Grid"]
        input_grid = CSV.File(input_grid_object["path"], stripwhitespace=true,  comment="#")

        # Required columns in grid file
        @assert :x in input_grid.names
        @assert :y in input_grid.names
        @assert :z in input_grid.names
        for row in input_grid
            mat = first([material for material in materials if material.id == row[:material]])
            push!(nodes, Nodes.Node(
                        Float64(row[:x]),
                        Float64(row[:y]),
                        Float64(row[:z]),
                        mat,
                        gridspacing
                    )
                 )
        end
    end

    # # Parse BCs
    # for bc in input["BC"]
        
    # end


    # Other (force planes, etc.)


    println("Finished reading input!")
end