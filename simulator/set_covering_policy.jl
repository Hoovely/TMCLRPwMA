function set_parameter(F::Vector{Int64}, G::Vector{Int64}, R::Vector{Int64}, TT::Matrix{Float64}, radius::Float64)
    
    I = vcat(G, R)
    J = F

    W = Dict()
    for i in I
        W[i] = [j for j in J if TT[i, j] <= radius]
    end

    return I, J, W
    
end

function set_covering_solve_location_problem(p::Int64, travel::Matrix{Float64}, F::Vector{Int64}, G::Vector{Int64}, R::Vector{Int64}, radius::Float64)
    
    I, J, W = set_parameter(F, G, R, travel, radius)

    println(W)

    # Initialize Cplex model
    model = Model(CPLEX.Optimizer)
    
    # Variable
    x = @variable(model, x[i in I], Bin)
    y = @variable(model, y[j in J], Bin)
    
    # Objective function
    @objective(model, Max, sum(x[i] for i in I))

    # Constraints
    # # Constraint 2
    # for i in I
    #     @constraint(model, x[i] == 1)
    # end

    # Constraint 3
    for i in I
        @constraint(model, sum(y[j] for j in W[i]) >= x[i])
    end

    # Constraint 4
    @constraint(model, sum(y[j] for j in J) == p)

    # Solve IP
    set_optimizer_attribute(model, "CPXPARAM_Threads", 15)
    # set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0)
    
    start_time = time()
    optimize!(model)
    solve_time = time() - start_time

    F_ = Vector{Int64}()
    for j in J 
        println(objective_value(model))
        if value.(y[j]) >= 0.9 
            push!(F_, j)
        end
    end

    return F_
end

function get_nearest_depot(node::Int64, F::Vector{Int64}, distance::Matrix{Float64})
    
    f′ = 0
    dis = 2147483647.0
    for f in F
        if dis > distance[f, node]
            dis = distance[f, node]
            f′ = f
        end
    end

    return f′
end

function set_covering_solve_assignment_problem(distance::Matrix{Float64}, F_::Vector{Int64}, G::Vector{Int64}, R::Vector{Int64})

    
    cluster_g = Dict{Int64, Vector}()
    for f in F_ cluster_g[f] = Vector{Int64}() end
    for g in G
        f′ = get_nearest_depot(g, F_, distance)
        push!(cluster_g[f′], g)
    end

    cluster_r = Dict{Int64, Vector}()
    for f in F_ cluster_r[f] = Vector{Int64}() end
    for r in R
        f′ = get_nearest_depot(r, F_, distance)
        push!(cluster_r[f′], r)
    end

    return cluster_g, cluster_r
end

function get_nearest_node(TT, service, now, next_set, visited)
    
    time = 2147483647.0
    node = 0
    for next in next_set
        if visited[next] == 1 continue end
        if time > TT[now, next] + service[next]
            time = TT[now, next] + service[next]
            node = next
        end
    end

    return node, time
end

function get_nearest_hospital(TT, demand, patient, hospital_set, hospital_capa)
    
    time = 2147483647.0
    hospital = 0
    for (idx, h) in enumerate(hospital_set)
        if demand > hospital_capa[idx] continue end
        if time > TT[patient, h]
            time = TT[patient, h]
            hospital = idx
        end
    end

    hospital_capa[hospital] -= demand

    return hospital, time, hospital_capa
end

function calculate_max_rescue_time(tours::Vector{Tour})
    
    rescue_time = 0
    for tour in tours
        rescue_time = max(rescue_time, tour.time)
    end

    return rescue_time
end

function set_covering_solve_routing_problem(node_size::Int64, F_::Vector{Int64}, H::Vector{Int64}, assigned_cluster_g::Dict{Int64, Vector}, assigned_cluster_r::Dict{Int64, Vector}, WG::Int64, WR::Int64,
    service::Vector{Int64}, patient::Vector{Int64}, demands::Vector{Int64}, distance::Matrix{Float64}, TT::Matrix{Float64}, capacity_v::Int64, capacity_man::Vector{Int64},
    EC_G::Int64, EC_R::Int64, VC_G::Int64, VC_R::Int64, FC::Int64)

    println(assigned_cluster_g)
    println(assigned_cluster_r)

    using_depot = Dict{Int64, Int64}()
    ne_tours, e_tours = Vector{Tour}(), Vector{Tour}()
    hospital_capa = copy(capacity_man)
    visited_g = Vector{Int64}(undef, node_size)
    visited_r = Vector{Int64}(undef, node_size)

    g_cost, r_cost = 0, 0
    for f_ in F_
        pre_node = 0
        patients_g = assigned_cluster_g[f_]
        patients_length = length(patients_g)

        flag, load, temp, capa, cnt = true, 0, 0, capacity_v, 1
        tour = Tour(Int[f_], 0.0, 0, Vector{Int}(), 0.0)
        while patients_length >= cnt
            if flag
                flag = false
                tour = Tour(Int[f_], 0.0, 0, Vector{Int}(), 0.0)

                node, t = get_nearest_node(TT, service, f_, patients_g, visited_g)    

                tour.operation_cost += distance[f_, node]
                capa -= demands[node]
                temp = demands[node]
            else
                node, t = get_nearest_node(TT, service, pre_node, patients_g, visited_g)

                if (node == 0) || (load + patient[node] > capa + temp - max(temp, demands[node]))
                    flag, load, temp, capa = true, 0, 0, capacity_v

                    tour.operation_cost += distance[node, f_] + VC_G + EC_G * tour.staff 
                    g_cost += tour.operation_cost
                    push!(tour.sequence, f_)
                    push!(tour.sequence, -1)
                    push!(ne_tours, tour)
                    using_depot[f_] = 1
                    continue
                end

                capa += temp
                temp = max(temp, demands[node])
                capa -= temp
                tour.operation_cost += distance[pre_node, node]
            end

            tour.time += t
            push!(tour.patient, patient[node])
            push!(tour.sequence, node)
            tour.staff = temp

            visited_g[node] = 1
            load += patient[node]
            pre_node = node

            cnt += 1
        end
        if flag == false
            tour.operation_cost += distance[tour.sequence[end], f_] + VC_G + EC_G * tour.staff 
            g_cost += tour.operation_cost

            push!(tour.sequence, f_)
            push!(tour.sequence, -1)
            push!(ne_tours, tour)
            using_depot[f_] = 1
        end


        patients_r = assigned_cluster_r[f_]
        patients_length = length(patients_r)

        cnt = 1
        while patients_length > cnt
            tour = Tour(Int[f_], 0.0, 0, Vector{Int}(), 0.0)

            r, t1 = get_nearest_node(TT, service, f_, patients_r, visited_r)
            push!(tour.patient, patient[r])
            push!(tour.sequence, r)
            tour.staff = demands[r]
            tour.time += t1
            visited_r[r] = 1

            h, t2, hospital_capa = get_nearest_hospital(TT, tour.staff, r, H, hospital_capa)
            push!(tour.sequence, node_size-length(H)+h)
            tour.time += t2
            tour.operation_cost = distance[f_, r] + distance[r, h] + VC_R + EC_R * tour.staff
            
            push!(e_tours, tour)
            cnt += 1
            r_cost += tour.operation_cost
            using_depot[f_] = 1
        end
    end

    println(ne_tours)
    println(e_tours)

    println(using_depot)

    total_cost = g_cost + r_cost + FC * length(collect(keys(using_depot)))
    obj = WG*calculate_max_rescue_time(ne_tours) + WR*calculate_max_rescue_time(e_tours)

    return ne_tours, e_tours, obj, total_cost

end