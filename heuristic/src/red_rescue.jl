"
    greedy_algorithm_for_red func
    Description: 
        - Calculating optimal time red patient to hospital
    Params:
        - Opening depot vector
        - Red patient vector
        - Hospital vector
        - Travel time matrix 
    Return:
        - Vector of red patient allocated hospital, Rescue completion time
"
function greedy_algorithm_for_red(using_depot_vec::Vector{Vector}, red_patient::Vector{Int64}, hospital::Vector{Int64}, travel_time::Matrix{Float64}, 
    distance::Matrix{Float64}, service::Vector{Int64}, patient::Vector{Int64}, demand::Vector{Int64}, capacity_room::Vector{Int64}, capacity_man::Vector{Int64}, VC_R::Int64, EC_R::Int64)

    obj_vec = Vector{Float64}()
    cost_vec = Vector{Float64}()
    trips_vec = Vector{Vector}()

    for using_depot in using_depot_vec

        obj = 0.0
        hospital_cnt = Dict()
        operation_cost = 0

        trips = Vector{Tour}(undef, length(red_patient))
        for (i, r) in enumerate(red_patient)
            tour = Vector{Int64}()
            rescue_completion_time = Inf
            total_cost = Inf

            for f in using_depot
                for h in hospital
                    total_time = travel_time[f, r] + travel_time[r, h] + service[r]
                    if rescue_completion_time > total_time
                        rescue_completion_time = total_time
                        tour = [f, r, h]
                        total_cost = distance[f, r] + distance[r, h] + VC_R + EC_R * demand[r]
                    end
                end
            end
            trips[i] = Tour(tour, rescue_completion_time, demand[r], Int[patient[r]], total_cost)
            obj = max(obj, trips[i].time)
            operation_cost += trips[i].operation_cost
        end

        push!(obj_vec, obj)
        push!(cost_vec, operation_cost)
        push!(trips_vec, trips)

    end

    return obj_vec, cost_vec, trips_vec
end

function greedy_algorithm_for_red_(using_depot::Vector{Int64}, red_patient::Vector{Int64}, hospital::Vector{Int64}, travel_time::Matrix{Float64}, 
    distance::Matrix{Float64}, service::Vector{Int64}, patient::Vector{Int64}, demand::Vector{Int64}, capacity_room::Vector{Int64}, capacity_man::Vector{Int64}, VC_R::Int64, EC_R::Int64)

    obj = 0.0
    hospital_cnt = Dict()
    operation_cost = 0

    trips = Vector{Tour}(undef, length(red_patient))
    for (i, r) in enumerate(red_patient)
        tour = Vector{Int64}()
        rescue_completion_time = Inf
        total_cost = Inf

        for f in using_depot
            for h in hospital
                total_time = travel_time[f, r] + travel_time[r, h] + service[r]
                if rescue_completion_time > total_time
                    rescue_completion_time = total_time
                    tour = [f, r, h]
                    total_cost = distance[f, r] + distance[r, h] + VC_R + EC_R * demand[r]
                end
            end
        end
        trips[i] = Tour(tour, rescue_completion_time, demand[r], Int[patient[r]], total_cost)
        obj = max(obj, trips[i].time)
        operation_cost += trips[i].operation_cost
    end

    return trips, operation_cost, obj
end


# ### Test Data ###

# # using Pkg
# # Pkg.add("Plots")
# # Pkg.add("LinearAlgebra")
# using Plots
# using LinearAlgebra  

# # Generate random (x, y) coordinates for 10 nodes
# node_coordinates = [(rand(0:10), rand(0:10)) for _ in 1:10]

# # Assign different nodes to depots, red patients, and hospitals
# depots = [5, 6]
# red_patients = [1, 2, 3, 4]
# hospitals = [7, 8, 9, 10]

# # Calculate Euclidean distances for the travel time matrix
# function calculate_distances(nodes)
#     num_nodes = length(nodes)
#     distances = zeros(num_nodes, num_nodes)
#     for i in 1:num_nodes
#         for j in 1:num_nodes
#             dx = nodes[i][1] - nodes[j][1]
#             dy = nodes[i][2] - nodes[j][2]
#             distances[i, j] = sqrt(dx^2 + dy^2)
#         end
#     end
#     distances
# end

# travel_time = calculate_distances(node_coordinates)


# # Run the algorithm
# chromosome, times = greedy_algorithm_for_red(depots, red_patients, hospitals, travel_time)
# println(node_coordinates)
# println(chromosome)
# print(times)

# function visualize_allocations(node_coordinates, depots, red_patients, hospitals, allocations)
#     plt = scatter([], [], labels="")  # Start with an empty plot to accumulate different scatter series

#     # Plot depots with labels
#     for d in depots
#         scatter!([node_coordinates[d][1]], [node_coordinates[d][2]], label="Depot $d", color=:yellow, shape=:square, markersize=8)
#         annotate!(node_coordinates[d][1] + 0.1, node_coordinates[d][2], text("Depot $d", 8, :left))
#     end

#     # Plot red patients with labels
#     for r in red_patients
#         scatter!([node_coordinates[r][1]], [node_coordinates[r][2]], label="Red Patient $r", color=:red, shape=:circle, markersize=8)
#         annotate!(node_coordinates[r][1] + 0.1, node_coordinates[r][2], text("Patient $r", 8, :left))
#     end

#     # Plot hospitals with labels
#     for h in hospitals
#         scatter!([node_coordinates[h][1]], [node_coordinates[h][2]], label="Hospital $h", color=:blue, shape=:diamond, markersize=8)
#         annotate!(node_coordinates[h][1] + 0.1, node_coordinates[h][2], text("Hospital $h", 8, :left))
#     end
    
#     # Draw paths
#     for (idx, r) in enumerate(red_patients)
#         depot, hospital = allocations[idx]
#         plot!([node_coordinates[depot][1], node_coordinates[r][1], node_coordinates[hospital][1]],
#               [node_coordinates[depot][2], node_coordinates[r][2], node_coordinates[hospital][2]],
#               label="", color=:green, line=(:dot, 2), arrow=true)
#     end
    
#     plot!(legend=:outertopright)
# end

# visualize_allocations(node_coordinates, depots, red_patients, hospitals, chromosome)
