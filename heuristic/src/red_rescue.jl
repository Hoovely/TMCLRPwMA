function set_parameter(TT::Matrix{Float64}, F::Vector{Int64}, R::Vector{Int64}, H::Vector{Int64}, service::Vector{Int64})

    cardinality_R = length(R)
    T = Matrix{Float64}(undef, length(R), length(F)*length(H))

    index_ = Dict()
    cnt = 1
    
    for h in H
        for f in F
            index_[cnt] = [f, h]
            
            for (idx, r) in enumerate(R)
                T[idx, cnt] = TT[f, r] + TT[r, h] + service[r]
            end
            cnt += 1
        end
    end
    
    I = [i for i in 1:cardinality_R]
    J = [j for j in 1:cnt-1]

    hospital_index_dict = Dict()
    lb = 1
    for h in H
        hospital_index_dict[h] = [f for f in J[lb:lb+length(F)-1]]
        lb += length(F)
    end

    return T, index_, I, J, hospital_index_dict
    
end

function resource_constraint_assignment_problem(TT::Matrix{Float64}, F::Vector{Int64}, R::Vector{Int64}, H::Vector{Int64}, service::Vector{Int64}, demand::Vector{Int64}, capacity_man::Vector{Int64}, slicing::Int64)

    T, index_, I, J, hospital_index_dict = set_parameter(TT, sort(F), R, H, service)

    # Initialize Cplex model
    model = Model(CPLEX.Optimizer)
    
    # Variable
    x = @variable(model, x[i in I, j in J], Bin)
    # x = @variable(model, x[i in I, j in J], lower_bound=0.0, upper_bound=1.0)
    eR = @variable(model, eR, lower_bound=0.0)

    # Objective function
    @objective(model, Min, eR)

    # Constraints
    # Constraint 2
    for i in I
        @constraint(model, sum(x[i, j] for j in J) == 1)
    end

    # # Constraint 3
    # for j in J
    #     @constraint(model, sum(x[i, j] for i in I) <= 1)
    # end

    # Constraint 4
    for h in H
        @constraint(model, sum(demand[i+slicing] * x[i, j] for j in hospital_index_dict[h] for i in I) <= capacity_man[h-slicing-length(R)])
    end

    # Constraint 5
    for i in I
        for j in J
            @constraint(model, T[i, j] * x[i, j] <= eR)
        end
    end

    # Solve IP
    set_optimizer_attribute(model, "CPXPARAM_Threads", 15)
    set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0)
    
    start_time = time()
    optimize!(model)
    solve_time = time() - start_time

    tour_vec = Vector{Vector}()
    for i in I 
        for j in J 
            if value.(x[i, j]) >= 0.9 
                # if value.(x[i, j]) != 1.0 println("Warning!") end
                f, h = index_[j][1], index_[j][2]
                push!(tour_vec, [f, i+slicing, h])
            end
        end
    end

    return objective_value(model), tour_vec
end

function LP_resource_constraint_assignment_problem(TT::Matrix{Float64}, F::Vector{Int64}, R::Vector{Int64}, H::Vector{Int64}, service::Vector{Int64}, demand::Vector{Int64}, capacity_man::Vector{Int64}, slicing::Int64)

    T, index_, I, J, hospital_index_dict = set_parameter(TT, sort(F), R, H, service)

    # Initialize Cplex model
    model = Model(CPLEX.Optimizer)
    
    # Variable
    x = @variable(model, x[i in I, j in J], lower_bound=0.0, upper_bound=1.0)
    # x = @variable(model, x[i in I, j in J], lower_bound=0.0, upper_bound=1.0)
    eR = @variable(model, eR, lower_bound=0.0)

    # Objective function
    @objective(model, Min, eR)

    # Constraints
    # Constraint 2
    for j in J
        @constraint(model, sum(x[i, j] for i in I) == 1)
    end

    # # Constraint 3
    # for j in J
    #     @constraint(model, sum(x[i, j] for i in I) <= 1)
    # end

    # Constraint 4
    for h in H
        @constraint(model, sum(demand[i+slicing] * x[i, j] for j in hospital_index_dict[h] for i in I) <= capacity_man[h-slicing-length(R)])
    end

    # Constraint 5
    for i in I
        for j in J
            @constraint(model, T[i, j] * x[i, j] <= eR)
        end
    end

    # Solve IP
    set_optimizer_attribute(model, "CPXPARAM_Threads", 15)
    set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0)
    
    start_time = time()
    optimize!(model)
    solve_time = time() - start_time

    tour_vec = Vector{Vector}()
    for i in I 
        for j in J 
            if value.(x[i, j]) >= 0.9 
                # if value.(x[i, j]) != 1.0 println("Warning!") end
                f, h = index_[j][1], index_[j][2]
                push!(tour_vec, [f, i+slicing, h])
            end
        end
    end

    return objective_value(model), tour_vec
end

function exact_assignment_algorithm_for_red(using_depot_vec::Vector{Vector}, F::Vector{Int64}, G::Vector{Int64}, R::Vector{Int64}, H::Vector{Int64}, TT::Matrix{Float64}, 
    distance::Matrix{Float64}, service::Vector{Int64}, patient::Vector{Int64}, demand::Vector{Int64}, capacity_room::Vector{Int64}, capacity_man::Vector{Int64}, VC_R::Int64, EC_R::Int64)

    obj_vec = Vector{Float64}()
    cost_vec = Vector{Float64}()
    trips_vec = Vector{Vector}()

    for using_depot in using_depot_vec

        operation_cost = 0
        trips = Vector{Tour}(undef, length(R))
        
        obj, tour_vec = resource_constraint_assignment_problem(TT, using_depot, R, H, service, demand, capacity_man, length(F)+length(G))

        for (idx, tour) in enumerate(tour_vec)
            f, r ,h = tour
            
            tour_time = TT[f, r] + TT[r, h] + service[r]
            total_cost = distance[f, r] + distance[r, h] + VC_R + EC_R * demand[r]

            trips[idx] = Tour(tour, tour_time, demand[r], Int[patient[r]], total_cost)

            operation_cost += trips[idx].operation_cost
        end

        push!(obj_vec, obj)
        push!(cost_vec, operation_cost)
        push!(trips_vec, trips)

    end

    return obj_vec, cost_vec, trips_vec
end


function exact_assignment_algorithm_for_red_(F::Vector{Int64}, R::Vector{Int64}, H::Vector{Int64}, TT::Matrix{Float64}, 
    distance::Matrix{Float64}, service::Vector{Int64}, patient::Vector{Int64}, demand::Vector{Int64}, capacity_room::Vector{Int64}, capacity_man::Vector{Int64}, VC_R::Int64, EC_R::Int64, slicing::Int64)

    operation_cost = 0
    trips = Vector{Tour}(undef, length(R))
    
    obj, tour_vec = resource_constraint_assignment_problem(TT, F, R, H, service, demand, capacity_man, slicing)

    for (idx, tour) in enumerate(tour_vec)
        f, r ,h = tour

        tour_time = TT[f, r] + TT[r, h] + service[r]
        total_cost = distance[f, r] + distance[r, h] + VC_R + EC_R * demand[r]

        trips[idx] = Tour(tour, tour_time, demand[r], Int[patient[r]], total_cost)

        operation_cost += trips[idx].operation_cost
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
