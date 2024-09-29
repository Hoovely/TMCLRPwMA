using Random
using StatsBase
using Combinatorics


function greedy_insertion_tour(T::Matrix{Float64}, depot::Int64, num_depot::Int64, t1::Vector{Int64})
    tour = Tour(Int[], 0.0)
    
    push!(tour.sequence, depot)

    current_node = depot
    while length(t1) > 0
        cost = Inf
        best_city_index = 0
        for (i, city) in enumerate(t1)
            if T[current_node, city+num_depot] < cost

                cost = T[current_node, city+num_depot]
                best_city_index = i
            end
        end
        push!(tour.sequence, t1[best_city_index]+num_depot)
        tour.cost += cost
        current_node = t1[best_city_index]+num_depot
        deleteat!(t1, best_city_index)
    end

    push!(tour.sequence, depot)
    tour.cost += T[current_node, depot]
    return tour
end

# Function to move the farthest median closer to the depot
# function move_farthest_median(medians, depot)
#     distances = [euclidean(median, depot) for median in medians]
#     farthest_median = argmax(distances)
#     medians[farthest_median] += 0.1 * (depot - medians[farthest_median])
#     return medians
# end


# function k_median(customers::Matrix{Float64}, depot::Vector{Float64}, k::Int)
#     data = customers
#     n = size(data)[1]
#     # Set the number of clusters you want to create and the depot location
#     # depot = [0.5, 0.5]

#     # Choose k initial cluster centers (medians)
#     medians = [data[i, :] for i in sample(1:n, k, replace=false)]

#     # Set a threshold for the difference between the maximum and minimum distances
#     # from the cluster medians to the depot
#     threshold = 0.1
#     assignments_ = Int[]
#     # Run the modified k-Median clustering algorithm
#     for i in 1:100
#         # Assign each node to the nearest median
#         assignments_ = [argmin([euclidean(data[p, :], m) for m in medians]) for p in 1:n]

#         # Compute the median of each cluster
#         for j in 1:k
#             cluster_points = data[findall(x -> x == j, assignments_), :]
#             if length(cluster_points) == 0
#                 medians[j] = depot
#             else
#                 medians[j] = median(cluster_points, dims=1)[1, :]
#             end
#         end

#         # Calculate the sum of distances from each cluster median to the depot
#         distances = [euclidean(median, depot) for median in medians]
#         max_distance = maximum(distances)
#         min_distance = minimum(distances)

#         # If the difference between the maximum and minimum distances is greater than
#         # the threshold, move the farthest median closer to the depot
#         if max_distance - min_distance > threshold
#             medians = move_farthest_median(medians, depot)
#         end
#     end
#     assignments_
# end

function calculate_center_point(t::Vector{Int64}, coor::Vector{Vector}, num_depot::Int64)
    cluster_x, cluster_y = 0, 0 
    for i in t
        cluster_x += coor[i+num_depot][1] 
        cluster_y += coor[i+num_depot][2]
    end 

    return [cluster_x/length(t), cluster_y/length(t)]
end

function calculate_nearest_depot(depot::Vector{Int64}, coor::Vector{Vector}, center_point::Vector{Float64})
    minimum_dis = Inf
    depot_node = -1

    for idx in depot
        dis = sqrt((center_point[1]-coor[idx][1])^2 + (center_point[2]-coor[idx][2])^2)
        
        if dis < minimum_dis
            minimum_dis = dis
            depot_node = idx
        end
    end


    return depot_node
end

function assign_patient(depot_vec::Vector{Vector}, Customers_::Matrix{Float64})

    Customers_ = transpose(Customers_)
    cluster_vec = Vector{Vector}()

    for open in depot_vec 
        number_open_depot = length(open)
        result = kmeans(Customers_, number_open_depot)
        assignments_ = copy(result.assignments)
        cluster_ = Vector{Vector{Int64}}(undef, number_open_depot)
        for i in 1:number_open_depot
            cluster_[i] = findall(x -> x == i, assignments_) # same clustering customer node index 
        end

        push!(cluster_vec, cluster_)
    end

    return cluster_vec
    
end

function allocate_depot(F::Vector{Int64}, assignment::Vector{Vector}, depot::Vector{Vector}, coordinate::Vector{Vector})

    tours_vec = Vector{Dict}()
    for (idx,comb) in enumerate(depot)
        number_open_depot = length(comb)
        tours_ = Dict{Int64, Vector{Int64}}()
        for i in 1:number_open_depot
            t1 = assignment[idx][i]
            center_point = calculate_center_point(t1, coordinate, length(F))
            depot_node = calculate_nearest_depot(comb, coordinate, center_point)

            if haskey(tours_, depot_node) == false tours_[depot_node] = Vector{Int64}() end
            for node in t1 push!(tours_[depot_node], node+length(F)) end
        end

        push!(tours_vec, tours_)
    end

    return tours_vec
end

# function initial_kmedian_solution(T::Matrix{Float64}, Customers_::Matrix{Float64}, depot::Vector{Int64}, opening_depot::Int64, coordinate::Vector{Vector},
#                                   distance::Matrix{Float64}, FC::Int64, EC_G::Int64, VC_G::Int64, r_cost::Float64, r_tour::Vector{Tour})
    
#     Customers_ = transpose(Customers_)

#     K = opening_depot 

#     if rand() < 1
#         result = kmeans(Customers_, K)
#         assignments_ = copy(result.assignments)
#     # else
#     #     assignments_ = k_median(Customers, depot, K)
#     end

#     tours = Vector{Tour}(undef, K)
#     genes = Int[]
#     obj = 0.0
#     for i in 1:K
#         t1 = findall(x -> x == i, assignments_) # same clustering customer node index 
#         center_point = calculate_center_point(t1, coordinate, length(depot))
#         depot_node = calculate_nearest_depot(depot, coordinate, center_point)
    
#         tour = greedy_insertion_tour(T, depot_node, length(depot), t1)
#         tours[i] = tour
#         genes = vcat(genes, tour.sequence)
#         genes = vcat(genes, [-1])

#         if tour.cost > obj
#             obj = tour.cost
#         end
#     end

#     chrm_cost = calculate_costs(trips, distance, FC, EC_G, VC_G) + r_cost
#     chrm = Chromosome(genes, obj, 0.0, chrm_cost, tours, r_tour)
#     return chrm
# end

# function initial_random_solution(T::Matrix{Float64}, K::Int, n_nodes::Int)
#     #     assignments_ = k_median(Customers, depot, K)

#     assignments_ = rand(1:K, n_nodes - K)
#     for i in 1:K
#         insert!(assignments_, rand(1:n_nodes-K+i), i)
#     end
#     tours = Tour[]
#     genes = Int[]
#     obj = 0.0
#     for i in 1:K
#         t1 = findall(x -> x == i, assignments_)
#         if length(t1) == 1
#             obj1 = T[1, t1[1]+1] + T[t1[1]+1, 1]
#             push!(tours, Tour(t1, obj1))
#             push!(genes, t1[1])
#             if obj1 > obj
#                 obj = obj1
#             end
#         else
#             t2 = copy(t1)
#             pushfirst!(t2, 0)
#             t2 = t2 .+ 1
#             TT = T[t2, t2]
#             tt1, obj1 = find_tsp_tour1(TT)
#             push!(tours, Tour(t1[tt1], obj1))
#             genes = vcat(genes, t1[tt1])
#             if obj1 > obj
#                 obj = obj1
#             end
#         end
#     end
#     chrm = Chromosome(genes, obj, 0.0, tours)
#     return chrm
# end

function find_tsp_tour(dist_mat::Matrix{Float64}, x_coor::Vector{Float64}, y_coor::Vector{Float64}, depot::Int64, tour::Vector{Int64}, method::String)

    if method == "Heuristic"
        result = TSP_Heuristic_solve_tsp(dist_mat, quality_factor=4)
        tsp_tour = result[1]
        tsp_tour[1], tsp_tour[end] = depot, depot
    end

    if method == "Hygese"
        result = Hygese_solve_tsp(x_coor, y_coor, AlgorithmParameters(timeLimit=0.1), verbose=false)
        tsp_tour = Vector{Int64}()
        push!(tsp_tour, depot)
        tsp_tour = vcat(tsp_tour, result.routes[1])
        push!(tsp_tour, depot)
    end

    if method == "LKH"
        # tsp_tour, _ = solve_tsp(x_coor, y_coor, dist="EUC_2D", log="off")
        # tour_length = length(tour)+2
    end

    tour_length = length(tour)+2

    for idx in 2:tour_length-1
        tsp_tour[idx] = tour[tsp_tour[idx]-1]
    end

    return tsp_tour[1:tour_length-1]
end

function make_matrix(depot::Int64, tour::Vector{Int64}, T::Matrix{Float64})
    
    node_size = length(tour)+1
    CT = zeros((node_size, node_size))
    for idx in 1:node_size
        for j in 1:node_size
            if idx == 1 || idx == node_size
                ## corner point ##
                if j == 1 || j == node_size
                    CT[idx, j] = 0.0 
                else
                    CT[idx, j] = T[depot, tour[j-1]]
                end
            else 
                if j == 1|| j == node_size
                    CT[idx, j] = T[depot, tour[idx-1]] ## side ##
                else
                    CT[idx, j] = T[tour[idx-1], tour[j-1]]
                end
            end
        end
    end

    return CT
end

function make_coor(depot::Int64, tour::Vector{Int64}, coordinate::Vector{Vector})
    x_coor = Vector{Float64}()
    y_coor = Vector{Float64}()

    push!(x_coor, Float64.(coordinate[depot][1]))
    push!(y_coor, Float64.(coordinate[depot][2]))
    for t in tour
        push!(x_coor, Float64.(coordinate[t][1]))
        push!(y_coor, Float64.(coordinate[t][2]))
    end
    # push!(x_coor, Float64.(coordinate[depot][1]))
    # push!(y_coor, Float64.(coordinate[depot][2]))

    return x_coor, y_coor
end



function find_mdtsp_tour(tours_vec::Vector{Dict}, distance::Matrix{Float64}, coordinate::Vector{Vector})
    
    mdtsp_tour_vec = Vector{Vector}()
    using_depot_vec = Vector{Vector}()

    for tours in tours_vec
        node_number_dict = Dict()
        mdtsp_tour = Vector{Vector{Int64}}()
        using_depot = Vector{Int64}()
        for key in keys(tours)
            dist_mat = make_matrix(key, tours[key], distance)
            x_coor, y_coor = make_coor(key, tours[key], coordinate)
            tsp_tour = find_tsp_tour(dist_mat, x_coor, y_coor, key, tours[key], "Heuristic")
            push!(mdtsp_tour, tsp_tour)
            push!(using_depot, key)
        end
        
        push!(mdtsp_tour_vec, mdtsp_tour)
        push!(using_depot_vec, using_depot)
    end
    return mdtsp_tour_vec, using_depot_vec
end


function open_depot(F::Vector{Int64}, mu::Int64)
    
    num_f = length(F)

    lb = 1
    ub = ceil(Int64, num_f*0.3)

    cnt = 0
    k_vec = Vector{Int64}()
    open_depot_vec = Vector{Vector}()
    for k in lb:ub
        push!(k_vec, k)
        for comb in collect(combinations(F, k))
            push!(open_depot_vec, comb)
            cnt = cnt + 1
        end
    end

    if cnt > mu
        open_depot_vec_ = Vector{Vector}()
        
        sample_size = ceil(Int64, mu/length(k_vec))

        e = 1
        l = 0
        for k in k_vec
            l += length(CoolLexCombinations(num_f, k))
            open_depot_vec_ = vcat(open_depot_vec_, sample(open_depot_vec[e:l], sample_size))
            e = l + 1
        end

        return open_depot_vec_, mu
    end

    return open_depot_vec, cnt
end

function check_depot(tours::Vector{Tour})

    depots = Vector{Int64}()
    for tour in tours
        push!(depots, tour.sequence[1])
    end

    depots = collect(Set(depots))

    return depots
end

# 현재는 node to node 변경만
# MDVRP 특성을 고려하여 다음과 같이 3가지 경우로 swap 가능
# 1. In cluster customer swap
# 2. Difference cluster customer swap
# 3. Depot to depot swap
# 이런 change를 통해 min-max를 개선하는 결과를 보이도록 휴리스틱 코드 짜야함


"
    change_initial func
    Description: 
        - Randomly change vehicles route sequence
    Params:
        - c::Vector{Int64} (Vehicle route for each depot, ex. depot1 customer1 customer2 customer3 customer4 customer5 ...)
        - n_nodes::Int64 (number of route length)
        - m::Int64 (number of vehicle)
    Return:
        - changed vehicles route sequence
"
function change_initial(c::Vector{Int64}, n_nodes::Int64, m::Int64)
    cc = copy(c)
    r = rand()
    if r < 0.4
        idx1 = rand(1:n_nodes)
        idx2 = rand(1:n_nodes)
        if idx1 > idx2
            temp = idx1
            idx1 = idx2
            idx2 = temp
        end
        cc[idx1:idx2] = reverse(cc[idx1:idx2])
    elseif r < 0.7
        if n_nodes <= 2 || n_nodes-1 < m return cc end
        idx = sort(sample(2:n_nodes-1, m - 1, replace=false))
        ccc = Vector{Vector{Int}}()
        k1 = 1
        for i in 1:m-1
            k2 = idx[i]
            push!(ccc, cc[k1:k2])
            k1 = k2 + 1
        end
        push!(ccc, cc[k1:n_nodes])
        shuffle!(ccc)
        cc = Int[]
        for t in ccc
            cc = vcat(cc, t)
        end
    else
        if n_nodes <= 2 return cc end
        idx = sample(1:n_nodes, rand(2:Int(ceil(n_nodes / 2))), replace=false)
        cc[idx] = shuffle(cc[idx])
    end
    return cc
end

function change_tour_to_gene(tours::Vector{Tour})
    genes = Vector{Int64}(undef, 0)
    for tour in tours
        genes = vcat(genes, tour.sequence)
        push!(genes, tour.sequence[1])
        push!(genes, -1)
    end

    return genes
end

function calculate_costs(trips::Vector{Tour}, r_tours::Vector{Tour}, distance::Matrix{Float64}, FC::Int64, EC_G::Int64, VC_G::Int64)

    total_cost = 0
    using_depot = Dict()
    for trip in trips
        length_trip = length(trip.sequence)
        using_depot[trip.sequence[1]] = 1
        trip.operation_cost += VC_G
        for (index, node) in enumerate(trip.sequence)
            if index == length_trip
                next_node = trip.sequence[1]
                trip.operation_cost += distance[node, next_node]
            else
                next_node = trip.sequence[index+1]
                trip.operation_cost += distance[node, next_node]
            end
        end
        trip.operation_cost += trip.staff * EC_G

        total_cost += trip.operation_cost
    end

    for r_tour in r_tours
        using_depot[r_tour.sequence[1]] = 1
    end

    opening_depot_cost = length(keys(using_depot))*FC
    
    total_cost += opening_depot_cost

    return total_cost
end

function create_random_chromosome(n_nodes::Int64)
    chromosome = shuffle!([i for i in 1:n_nodes])
    chromosome
end
function create_random_chromosome2(T::Matrix{Float64}, F::Vector{Int64}, G::Vector{Int64}, n_nodes::Int64, m::Int)
    customers = shuffle([i for i in G])
    tours = Tour[]
    for i in 1:m
        push!(tours, Tour(Int[rand(F)], 0.0))
    end

    while length(customers) > 0
        city = pop!(customers)
        put_city_in_tour(tours, city, T, n_nodes)
    end

    S = Int[]
    obj = 0.0
    for tour in tours
        S = vcat(S, tour.sequence)
    end
    return S
end

function create_random_chromosome3(T::Matrix{Float64}, n_nodes::Int64)
    Nodes = [i for i in 1:n_nodes]
    tours = Tour[]
    a = copy(T[1, 2:n_nodes+1])
    b = sortperm(a)

    for i in 1:m
        push!(tours, Tour([b[i]], 0.0))
    end

    tour_indices = [i for i in 1:m]
    while length(Nodes) > 0
        r = tour_indices[rand(1:length(tour_indices))]
        tour = tours[r]
        last_city = tour.sequence[length(tour.sequence)]
        a = copy(T[last_city+1, Nodes.+1])
        new_city = Nodes[argmin(a)]
        push!(tour.sequence, new_city)
        deleteat!(tour_indices, findfirst(x -> x == r, tour_indices))
        deleteat!(Nodes, findfirst(x -> x == new_city, Nodes))
        if length(tour_indices) == 0
            tour_indices = [i for i in 1:m]
        end
    end

    S = Int[]
    for tour in tours
        S = vcat(S, tour.sequence)
    end
    return S
end

function generate_initial_population(TT::Matrix{Float64}, mu::Int64, sigma::Int64, num_depot_vec::Int64, tsp_tours_vec::Vector{Vector}, r_obj::Vector{Float64}, r_cost::Vector{Float64}, r_tour::Vector{Vector}, Customers_::Matrix{Float64}, 
    F::Vector{Int64}, G::Vector{Int64}, R::Vector{Int64}, H::Vector{Int64}, Return_Depot::Vector{Int64}, distance::Matrix{Float64}, service::Vector{Int64}, patient::Vector{Int64},
    demand::Vector{Int64}, capacity_vg::Int64, capacity_vr::Int64, FC::Int64, EC_G::Int64, EC_R::Int64, VC_G::Int64, VC_R::Int64, WG::Int64, WR::Int64)

    population = Vector{Chromosome}(undef, sigma)

    println("Start Initialization")

    for (idx, tsp_tours) in enumerate(tsp_tours_vec)
        obj, trips = SPLIT(TT, tsp_tours, service, patient, demand, capacity_vg)
        chrm_cost = calculate_costs(trips, r_tour[idx], distance, FC, EC_G, VC_G)
        best_mdvrp_tours = change_tour_to_gene(trips)
        population[idx] = Chromosome(best_mdvrp_tours, WG*obj+WR*r_obj[idx], obj, 0.0, chrm_cost+r_cost[idx], check_depot(trips), trips, r_tour[idx], r_cost[idx], r_obj[idx])
    end

    sort!(population[1:num_depot_vec], by=x -> x.weighted_time)

    for i in 1+num_depot_vec:sigma
        j = i-num_depot_vec
        if rand() < 1
            changed_tsp_tours = Vector{Vector{Int64}}()
            for tour in population[j].tours
                depot, S = tour.sequence[1], tour.sequence[2:end]
                K = Get_feasible_K(S, capacity_vg, patient, demand)
                push!(changed_tsp_tours, vcat(depot, change_initial(S, length(S), K)))
            end
            
            obj, trips = SPLIT(TT, changed_tsp_tours, service, patient, demand, capacity_vg)
            chrm_cost = calculate_costs(trips, population[j].reds_path, distance, FC, EC_G, VC_G)
            S = change_tour_to_gene(trips)
            population[i] = Chromosome(S, WG*obj+WR*population[j].r_obj, obj, 0.0, chrm_cost+population[j].r_cost, population[j].depots, trips, population[j].reds_path, population[j].r_cost, population[j].r_obj)
        else  
            # population[i] = initial_kmedian_solution(TT, Customers_, depot_vec, number_open_depot, coordinate, distance, FC, EC_G, VC_G, r_cost, r_tour)
        end
    end

    sort!(population, by=x -> x.weighted_time)

    return population, population[1].weighted_time
    
end



"
    diversify! func
    Description: 
        - Diversify for population, randomly change to MDVRP route by using some functions (change_initial, create_random_chromosome2, initial_kmedian_solution)
    Params:
        - population::Vector{Chromosome} (population set)
        - TT::Matrix{Float64} (distance matrix about all of nodes)
        - K::Int (number of vehicles) 
        - mu::Int (number of populations) 
        - tsp_tours::Vector{Vector{Int}} (tsp tours for each depots) 
        - customers::Matrix{Float64} () 
        - depot::Vector{Float64} ()
    Return:
        - changed vehicles route sequence
"
function diversify!(population::Vector{Chromosome}, TT::Matrix{Float64}, K::Int, mu::Int, tsp_tours::Vector{Vector{Int}}, 
    customers::Matrix{Float64}, depot::Vector{Int64}, opening_depot::Int64, coordinate::Vector{Vector{Int64}}, F::Vector{Int64}, G::Vector{Int64}, num::Int)
    # println("Start diversify!")

    n_tours = length(tsp_tours)
    n_nodes = size(TT)[1] - 2
    n_best = Int(round(0.15 * mu)) # number of save best chromosome?
    for i = n_best+1:length(population) # length(population) == mu
        if rand() < 0.8
            SS = Vector{Vector{Int64}}()
            if rand() < 0.7
                r = rand(1:n_tours) # select to changing depot
                tour = tsp_tours[r]
                scaled_K = Int64(round(K*length(tour)/length(G)))
                push!(SS, vcat(tour[1], change_initial(tour[2:end], length(tour[2:end]), scaled_K)))
            else
                push!(SS, create_random_chromosome2(TT, F, G, n_nodes, K))
            end

            obj, trips = SPLIT(TT, K, SS, length(G))
            S = change_tour_to_gene(trips)
            chrm = Chromosome(S, obj, 0.0, trips)
        else
            chrm = initial_kmedian_solution(TT, customers, depot, opening_depot, coordinate)
        end
        population[i] = chrm
    end
    sort!(population, by=x -> x.fitness)
end


function find_tour_neighbors(tours::Vector{Tour}, Customers::Matrix{Float64}, depot::Vector{Float64}, m::Int)
    means = [mean(Customers[t1, :], dims=1)[1, :] for t1 in [tours[i].sequence for i in 1:m]]
    for (j, v) in enumerate(means)
        means[j] = (m * v + depot) / (m + 1)
    end
    distances = ones(m, m) * -1
    for i in 1:m-1
        for j = i+1:m
            distances[i, j] = euclidean(means[i], means[j])
            distances[j, i] = distances[i, j]
        end
    end
    return [sortperm(distances[i, :])[2:min(m, 10)] for i in 1:m]
end

function find_tour_neighbors(tours::Vector{Tour}, T::Matrix{Float64}, m::Int)

    distances = ones(m, m) * -1
    # distance check, tour[i] and tour[j]
    for i in 1:m-1
        for j = i+1:m
            distances[i, j] = sum(T[tours[i].sequence[2:end], tours[j].sequence[2:end]]) # 
            distances[j, i] = distances[i, j]
        end
    end

    # return to neighborhood tours min(m, 5)
    return [sortperm(distances[i, :])[2:min(m, 5)] for i in 1:m]
end


function enrich_the_chromosome!(chrm::Chromosome, T::Matrix{Float64}, customers::Matrix{Float64}, depot::Vector{Float64}, n_nodes::Int)

    m = length(chrm.tours)
    temp = deepcopy(chrm)
    max_tour_index = argmax([chrm.tours[i].cost for i in 1:length(chrm.tours)])
    max_tour_length = chrm.fitness
    #     tour_neighbors = find_tour_neighbors(Chrm.tours, Customers, depot, m)
    tour_neighbors = find_tour_neighbors(chrm.tours, T, m)
    improved = true
    count = 0

    while improved && count < 100
        count += 1
        improved = false
        for r1  in 1:m
            for r2 in tour_neighbors[r1]
                if r2 != max_tour_index

                    tour1 = chrm.tours[r1].sequence
                    tour2 = chrm.tours[r2].sequence

                    k1 = 0
                    became_max = false
                    count2 = 0
                    while k1 < length(tour1) && !became_max
                        
                        cost1 = chrm.tours[r1].cost
                        cost2 = chrm.tours[r2].cost
                        k1 += 1
                        city1 = tour1[k1]
                        for k2  in 1:length(tour2)+1
                            if chrm.tours[r1].cost > chrm.tours[r2].cost
                                new_cost2 = calculate_new_cost_add_one(tour2, cost2, city1, k2, T, n_nodes)
                                new_cost1 = calculate_new_cost_remove_one(tour1, cost1, k1, T, n_nodes)
                                do_it = false
                                if r1 == max_tour_index
                                    if new_cost2 < max_tour_length
                                        do_it = true
                                    end
                                else
                                    if (new_cost2 - cost2) < (cost1 - new_cost1) && new_cost2 < max_tour_length
                                        do_it = true
                                    end
                                end
                                if do_it

                                    chrm.tours[r1].cost = new_cost1
                                    chrm.tours[r2].cost = new_cost2
                                    insert!(tour2, k2, city1)
                                    deleteat!(tour1, k1)
                                    if new_cost2 > max_tour_length
                                        max_tour_length = new_cost2
                                        max_tour_index = r2
                                        became_max = true
                                    end
                                    k1 -= 1
                                    improved = true
                                    break
                                end
                            end
                        end
                    end
                end
            end
            
        end
    end
    for tour in chrm.tours
        two_opt_on_route(tour, T, n_nodes)
    end
    chrm.genes = Int[]
    chrm.fitness = maximum([chrm.tours[i].cost for i in 1:length(chrm.tours)])
    for tour in chrm.tours
        chrm.genes = vcat(chrm.genes, tour.sequence)
    end
end

function find_closeness(TT::Matrix{Float64}, num_depot::Int64, num_patient::Int64, h::Float64)
    # find closeness nodes matrix for each node 

    num = Int(ceil(h * num_patient))
    closenessT = fill(false, num_patient + num_depot, num_patient + num_depot)

    for i = num_depot+1:num_depot+num_patient
        a = copy(TT[i, num_depot+1:num_depot+num_patient])
        b = sortperm(a)
        closenessT[i, b[2:num+1]] .= true
    end
    return closenessT
end