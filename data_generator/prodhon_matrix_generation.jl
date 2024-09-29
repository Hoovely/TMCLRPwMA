function prodhon_matrix(coordinates::Vector{Vector}, depot_num::Int64)
    
    num_of_nodes = length(coordinates)
    dist_matrix = zeros(Float64, num_of_nodes+depot_num, num_of_nodes+depot_num)
    
    for i in 1:num_of_nodes
        temp = coordinates
        for j in 1:num_of_nodes
            distance = sqrt((temp[i][1]-temp[j][1])^2+(temp[i][2]-temp[j][2])^2)
            dist_matrix[i, j] = distance
            if i <= depot_num
                dist_matrix[i+num_of_nodes, j] = distance
            end
            if j <= depot_num
                dist_matrix[i, j+num_of_nodes] = distance
            end
        end
    end

    return dist_matrix
end