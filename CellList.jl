struct CellList
    data::Array{Vector{Main.PD.Node}, 3}
    division_size::Float64
    data_min::Vector{Float64}
    data_max::Vector{Float64}
end

function create_cell_list(nodes::Vector{Main.PD.Node}, radius::Float64)
    # Can be optimized
    positions = [node.position + node.displacement for node in nodes]
    dim_max = copy(positions[1])
    dim_min = copy(positions[1])
    for position in positions
        if position[1] > dim_max[1]
            dim_max[1] = position[1]
        end
        if position[2] > dim_max[2]
            dim_max[2] = position[2]
        end
        if position[3] > dim_max[3]
            dim_max[3] = position[3]
        end
        if position[1] < dim_min[1]
            dim_min[1] = position[1]
        end
        if position[2] < dim_min[2]
            dim_min[2] = position[2]
        end
        if position[3] < dim_min[3]
            dim_min[3] = position[3]
        end
    end
    shape = (dim_max - dim_min) / radius
    size1,size2,size3 = Int64(ceil(shape[1]))+1, Int64(ceil(shape[2]))+1,Int64(ceil(shape[3]))+1

    data = Array{Vector{Main.PD.Node}}(undef,size1,size2,size3);
    for i in eachindex(data)
        data[i] = Vector{Main.PD.Node}()
    end


    for node in nodes
        # Insert node
        pos = node.position + node.displacement
        # index = []
        push!(data[Int64(ceil((pos[1] - dim_min[1]) / radius))+1,
                    Int64(ceil((pos[2] - dim_min[2]) / radius))+1,
                    Int64(ceil((pos[3] - dim_min[3]) / radius))+1], node)
    end

        
    println("Created cell list: ", size(data))
    return CellList(data, radius, dim_min, dim_max)
end

function create_cell_list_new(nodes::Vector{Main.PD.Node}, radius::Float64)
    # Can be optimized
    positions = [node.position + node.displacement for node in nodes]
    dim_max = copy(positions[1])
    dim_min = copy(positions[1])
    for position in positions
        if position[1] > dim_max[1]
            dim_max[1] = position[1]
        end
        if position[2] > dim_max[2]
            dim_max[2] = position[2]
        end
        if position[3] > dim_max[3]
            dim_max[3] = position[3]
        end
        if position[1] < dim_min[1]
            dim_min[1] = position[1]
        end
        if position[2] < dim_min[2]
            dim_min[2] = position[2]
        end
        if position[3] < dim_min[3]
            dim_min[3] = position[3]
        end
    end
    shape = (dim_max - dim_min) / radius
    size1,size2,size3 = Int64(ceil(shape[1]))+1, Int64(ceil(shape[2]))+1,Int64(ceil(shape[3]))+1

    a = Vector{Array{Vector{Main.PD.Node},3}}()
    for i in 1:Threads.nthreads()
        b = Array{Vector{Main.PD.Node}}(undef,size1,size2,size3)
        for j in eachindex(b)
            b[j] = Vector{Main.PD.Node}()
        end
        push!(a, b)
    end

    Threads.@threads for node in nodes
        # Insert node
        pos = node.position + node.displacement
        index1, index2, index3 = Int64(ceil((pos[1] - dim_min[1]) / radius))+1, Int64(ceil((pos[2] - dim_min[2]) / radius))+1, Int64(ceil((pos[3] - dim_min[3]) / radius))+1
        push!(a[Threads.threadid()][index1, index2, index3], node)
    end
    return CellList(vcat(a...), radius, dim_min, dim_max)
end

function sample_cell_list(cell_list::CellList, node::Main.PD.Node, radius::Float64)
    pos = node.position + node.displacement
    index = [
        Int64(ceil((pos[1] - cell_list.data_min[1]) / radius))+1,
        Int64(ceil((pos[2] - cell_list.data_min[2]) / radius))+1,
        Int64(ceil((pos[3] - cell_list.data_min[3]) / radius))+1
        ]   
    shape = size(cell_list.data)
    result = Vector{Main.PD.Node}()

    for i in -1:1
        if index[1]+i > 0 && index[1]+i <= shape[1]
            for j in -1:1
                if index[2]+j > 0 && index[2]+j <= shape[2]
                    for k in -1:1
                        if index[3]+k > 0 && index[3]+k <= shape[3]
                            append!(result, cell_list.data[index[1]+i, index[2]+j, index[3]+k])
                        end
                    end
                end
            end
        end
    end
    return [n for n in result if norm(n.position+n.displacement-node.position-node.displacement) < radius && n != node]
end


# nodes = Vector{Main.PD.Node}()
# mat = PD.Material(1, 1., 0.1, 75)

# push!(nodes, Main.PD.Node(-1.5,-1.1,-1.3,mat, 1.))
# push!(nodes, Main.PD.Node(1.4,1.6,1.8,mat, 1.))

# @time create_cell_list(nodes, 0.7);
# @time create_cell_list_new(nodes, 0.7);

# Threads.nthreads()




# a = Main.PD.Node(1.,0.,0.,mat, 1.)
# sample_cell_list(cell_list, a, 1.314)






# push!(nodes, Main.PD.Node(0.,0.,0.,mat, 1.))
# push!(nodes, Main.PD.Node(1.,0.,0.,mat, 1.))
# push!(nodes, Main.PD.Node(2.,0.,0.,mat, 1.))
# push!(nodes, Main.PD.Node(3.,0.,0.,mat, 1.))
# create_cell_list(nodes)

