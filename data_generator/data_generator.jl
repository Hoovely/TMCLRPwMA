using Random

include("augerat_matrix_generation_cut.jl")
include("load_data_augerat.jl")
include("duplicate_data_augerat.jl")


Random.seed!(0)


function get_data(original_data::Dict{Any, Any}, problem::String, augerat_data_::Dict{Any, Any}, cutting_num::Int64, tmc_num::Int64)
    # Matrix Generation
    data_tot = create_node_data(problem, augerat_data_, cutting_num) # cutting_num parameter 값에 따라 node cutting. node의 coordinates, qty(demand) return
    mat = original_data[problem]
    M, M_data = dup_data(data_tot, mat, cutting_num, tmc_num)

    return M, M_data
end

function generate_hospital(coor::Vector{Vector}, hosptital_num::Int64)

    min_x, min_y = 2147483647, 2147483647
    max_x, max_y = 1, 1
    for (x, y) in coor
        min_x, min_y = min(min_x, x), min(min_y, y)
        max_x, max_y = max(max_x, x), max(max_y, y)
    end
    
    hos_coordi = Vector{Vector}()
    for _ in 1:hosptital_num
        push!(hos_coordi, [ rand(min_x:max_x), rand(min_y:max_y) ])
    end

    return hos_coordi
end

function generate_node(idx_list::Vector{Int64}, number::Int64)
    return [idx_list[i] for i in eachindex(idx_list) if i <= number]
end

function generate_staff(idx_list::Vector{Int64}, facility_list::Vector{Int64}, patient_number::Int64)
    for i in eachindex(idx_list)
        if i in facility_list 
            idx_list[i] = trunc(Int64, patient_number/length(facility_list))*5 
        end
    end
    return idx_list
end

function generate_patient(idx_list::Vector{Int64}, G::Vector{Int64}, R::Vector{Int64})
    for i in eachindex(idx_list)
        if i in G
            idx_list[i] = rand(1:3)
        end

        if i in R
            idx_list[i] = 1
        end
    end
    return idx_list
end 

function generate_demand(idx_list::Vector{Int64}, patient_list::Vector{Int64})
    for i in eachindex(idx_list)
        if i in patient_list
            idx_list[i] = rand(1:3)
        end
    end
    return idx_list
end 

function difference(lst1::Vector{Int64}, lst2::Vector{Int64})
    return sort( collect( setdiff(Set(lst1), Set(lst2)) ) )
end

function generate_service(N::Vector{Int64}, facility::Vector{Int64}, green::Vector{Int64}, red::Vector{Int64}, hospital::Vector{Int64})
    vec = Vector{Int64}()
    for i in N
        if i in facility 
            push!(vec, 0) 
        elseif i in green
            push!(vec, rand(5:35))
        elseif i in red
            push!(vec, rand(2:15))
        elseif i in hospital
            push!(vec, 0)
        end
    end
    return vec
end

function generate_travel(distance_matrix::Matrix{Float64}, velocity::Int64, scale::Int64)
    n = size(distance_matrix)[1]
    time_matrix = zeros(n, n)
    for i in 1:n
        for j in 1:n
            time_matrix[i, j] = distance_matrix[i, j]*scale/velocity 
        end
    end
    return time_matrix
end

function generate_h_room(patient::Vector{Int64}, R::Vector{Int64}, H::Vector{Int64})
    
    red_sum = 0
    for (idx,p) in enumerate(patient)
        if idx in R 
            red_sum += p 
        end
    end

    capa = ceil(Int64, red_sum/length(H))
    
    vec = Vector{Int64}()
    for _ in H
        push!(vec, capa)
    end

    return vec
end

function generate_h_man(demand::Vector{Int64}, R::Vector{Int64}, H::Vector{Int64})
    
    red_sum = 0
    for (idx,d) in enumerate(demand)
        if idx in R 
            red_sum += d
        end
    end

    capa = ceil(Int64, red_sum/length(H))
    
    vec = Vector{Int64}()
    for _ in H
        push!(vec, capa)
    end

    return vec
end