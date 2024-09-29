function new_build_vehicle(new_chrm::Chromosome, TT::Matrix{Float64}, distance::Matrix{Float64}, service::Vector{Int64}, patient::Vector{Int64}, 
                           demand::Vector{Int64}, FC::Int64, EC_G::Int64, VC_G::Int64, epsilon::Int64, r_tours::Vector{Tour}, r_cost::Float64)

    single_cnt = 1
    while epsilon >= new_chrm.total_cost
        if single_cnt > length(new_chrm.tours) break end

        g_obj, index = sort([(tour.time, idx) for (idx,tour) in enumerate(new_chrm.tours)], rev=true)[single_cnt]
        f, S = new_chrm.tours[index].sequence[1], copy(new_chrm.tours[index].sequence[2:end])

        if length(S) == 1
            single_cnt += 1
            continue
        end

        cut_point = round(Int64, length(S)/2)

        new_tour = Vector{Tour}()
        time1 = 0.0
        tour1 = Tour(vcat(f, S[1:cut_point]), time1, 0, Int[], 0.0)
        for (idx, node) in enumerate(tour1.sequence)
            if node == tour1.sequence[end]
                next = tour1.sequence[1]
            else
                next = tour1.sequence[idx+1]
            end

            time1 += TT[node, next] + service[node]

            if idx > 1
                tour1.staff = max(tour1.staff, demand[node])
                push!(tour1.patient, patient[node])
            end
        end
        tour1.time = time1
        push!(new_tour, tour1)

        time2 = 0.0
        tour2 = Tour(vcat(f, S[cut_point+1:end]), time2, 0, Int[], 0.0)
        for (idx, node) in enumerate(tour2.sequence)
            if node == tour2.sequence[end]
                next = tour2.sequence[1]
            else
                next = tour2.sequence[idx+1]
            end
            
            time2 += TT[node, next] + service[node]

            if idx > 1
                tour2.staff = max(tour2.staff, demand[node])
                push!(tour2.patient, patient[node])
            end
        end
        tour2.time = time2
        push!(new_tour, tour2)

        for (idx, tour) in enumerate(new_chrm.tours)
            if idx == index continue end
            tour.operation_cost = 0.0
            push!(new_tour, tour)
        end

        updated_cost = calculate_costs(new_tour, r_tours, distance, FC, EC_G, VC_G)+r_cost

        if updated_cost >= epsilon break end

        new_chrm.tours = new_tour

        single_cnt = 1
    end
    
    # if new_chrm.total_cost == 0.0
    #     println("In Rebalancing")
    #     println(new_chrm)
    #     wait()
    # end

    return new_chrm
end


function evaluate_f′_r(TT::Matrix{Float64}, service::Vector{Int64}, r_tours::Vector{Tour}, f′::Int64)

    improvement = 0.0
    for r_tour in r_tours
        f, r, h = r_tour.sequence

        diff = r_tour.time - (TT[f′, r] + service[r] + TT[r, h])
        if diff > 0
            improvement += diff
        end
    end

    return improvement
end


function evaluate_f′_g(TT::Matrix{Float64}, depots_::Matrix{Float64}, F_::Vector{Int64}, f′::Int64)

    if rand() > 0.5
        return mean([TT[f, f′] for f in F_])
    else
        x, y = mean([depots_[f, 1] for f in F_]), mean([depots_[f, 2] for f in F_])
        x′, y′ = depots_[f′, 1], depots_[f′, 2]
        return euclidean([x,y], [x′, y′])
    end
end


function new_open_depot(chrm_::Chromosome, TT::Matrix{Float64}, service::Vector{Int64}, depots_::Matrix{Float64}, WG::Int64, WR::Int64, F_′::Vector{Int64}, F′::Vector{Int64})
  
    score = zeros(length(F′))
    for (idx,f′) in enumerate(F′)
        # score in emergency
        score_r = evaluate_f′_r(TT, service, chrm_.reds_path, f′)

        # score in non-emergency
        score_g = evaluate_f′_g(TT, depots_, F_′, f′)

        # println("red score is $(score_r)")
        # println("green score is $(score_g)")

        score[idx] = WR*score_r/(score_r+score_g) + WG*score_g/(score_r+score_g)
    end
    
    f_′ = F′[argmax(score)]
    push!(F_′, f_′)

    return F_′
end

function new_clustering_routing_split(F::Vector{Int64}, depot::Vector{Int64}, Customers_::Matrix{Float64}, coordinate::Vector{Vector}, distance::Matrix{Float64},
                                      TT::Matrix{Float64}, service::Vector{Int64}, patient::Vector{Int64}, demand::Vector{Int64}, capacity_vg::Int64)

    # Clustering
    Customers_ = transpose(Customers_)
    number_opening_depot = length(depot)
    result = kmeans(Customers_, number_opening_depot)
    assignments_ = copy(result.assignments)
    cluster_ = Vector{Vector{Int64}}(undef, number_opening_depot)
    for i in 1:number_opening_depot
        cluster_[i] = findall(x -> x == i, assignments_) # same clustering customer node index 
    end

    # Allocation
    tours = Dict{Int64, Vector{Int64}}()
    for i in 1:number_opening_depot
        t1 = cluster_[i]
        center_point = calculate_center_point(t1, coordinate, length(F))
        depot_node = calculate_nearest_depot(depot, coordinate, center_point)

        if haskey(tours, depot_node) == false tours[depot_node] = Vector{Int64}() end
        for node in t1 push!(tours[depot_node], node+length(F)) end
    end

    # Routing
    mdtsp_tour = Vector{Vector{Int64}}()
    for key in keys(tours)
        dist_mat = make_matrix(key, tours[key], distance)
        x_coor, y_coor = make_coor(key, tours[key], coordinate)
        tsp_tour = find_tsp_tour(dist_mat, x_coor, y_coor, key, tours[key], "Heuristic")
        push!(mdtsp_tour, tsp_tour)
    end

    # SPLIT
    obj, trips = SPLIT(TT, mdtsp_tour, service, patient, demand, capacity_vg)

    return obj, trips
end

function rebalancing(population::Vector{Chromosome}, TT::Matrix{Float64}, customers_::Matrix{Float64}, depots_::Matrix{Float64}, F::Vector{Int64}, G::Vector{Int64},
                      R::Vector{Int64}, H::Vector{Int64}, distance::Matrix{Float64}, service::Vector{Int64}, patient::Vector{Int64}, demand::Vector{Int64}, 
                      capacity_vg::Int64, capacity_vr::Int64, capacity_room::Vector{Int64}, capacity_man::Vector{Int64}, FC::Int64, EC_G::Int64, EC_R::Int64, VC_G::Int64, VC_R::Int64, WG::Int64, WR::Int64,
                      coordinate::Vector{Vector}, epsilon::Int64)

    for (idx, chrm) in enumerate(population)
        chrm_ = deepcopy(chrm)
        F_ = chrm_.depots
        operation_cost = chrm_.total_cost
        tours_ = chrm_.tours
        r_tours, r_obj, r_cost = chrm_.reds_path, chrm_.r_obj, chrm_.r_cost
        remain = epsilon-operation_cost

        if remain <= 0 continue end

        add_num_g = floor(Int64, length(G)/(length(F_)+1))
        K′ = ceil(Int64, add_num_g/5)
        exp_cost = (FC + K′ * VC_G) * 1.2

        if remain < exp_cost || length(F_) == length(F)
            max_route_length = 0
            for tour in tours_
                max_route_length = max(max_route_length, length(tour.sequence))
            end
            if max_route_length == 2 continue end

            # Add vehicle
            chrm = Chromosome(Int[], 0.0, 0.0, 0.0, 0.0, F_, tours_, r_tours, r_cost, r_obj)
            new_chrm = new_build_vehicle(chrm_, TT, distance, service, patient, demand, FC, EC_G, VC_G, epsilon, r_tours, r_cost)

            obj = 0.0
            for tour in new_chrm.tours 
                obj = max(obj, tour.time) 
            end 

            new_chrm.genes = change_tour_to_gene(new_chrm.tours)
            new_chrm.weighted_time = WG*obj+WR*r_obj
            new_chrm.fitness = obj

            population[idx] = new_chrm
        else
            # Build depot
            F′ = collect(setdiff(Set(F), Set(F_)))
            F_′ = new_open_depot(chrm_, TT, service, depots_, WG, WR,  copy(F_), F′)
            obj′, trips = new_clustering_routing_split(F, F_′, customers_, coordinate, distance, TT, service, patient, demand, capacity_vg)
            tours_′ = change_tour_to_gene(trips)

            # Calculate total cost
            r_tours′, r_obj′, r_cost′ = greedy_algorithm_for_red_(F_′, R, H, TT, distance, service, patient, demand, capacity_room, capacity_man, VC_R, EC_R)
            operation_cost_ = calculate_costs(trips, r_tours, distance, FC, EC_G, VC_G)+r_cost

            if operation_cost_ > epsilon
                # println("============add open cost============")
                # println("Over the remain budget!!! $(operation_cost_) > $(remain)")
                continue
            end

            population[idx] = Chromosome(tours_′, WG*obj′+WR*r_obj, obj′, 0.0, operation_cost_, F_′, trips, r_tours′, r_obj′, r_cost′)
        end
    end

    return population
end