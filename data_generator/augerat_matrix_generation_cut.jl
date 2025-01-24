function augerat_matrix(augerat_data_::Dict{Any, Any}, name::String)
    
    dist_mat_dic = Dict()
    num_of_nodes = length(augerat_data_[name]["coordinates"])
    dist_matrix = Matrix{Float64}(undef, num_of_nodes, num_of_nodes)
    
    for i in 1:num_of_nodes
        # push!(dist_matrix, [])
        
        temp = augerat_data_[name]["coordinates"]
        for j in 1:num_of_nodes
            distance = sqrt((temp[string(i)][1]-temp[string(j)][1])^2 
                            +(temp[string(i)][2]-temp[string(j)][2])^2)
            dist_matrix[i, j] = distance
            # push!(dist_matrix[i], distance)
        end
        
    end
    dist_mat_dic[name]=dist_matrix

    return dist_mat_dic
end