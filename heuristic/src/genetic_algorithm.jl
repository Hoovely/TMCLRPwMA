mutable struct Tour
    sequence::Vector{Int64}
    time::Float64
    staff::Int64
    patient::Vector{Int64}
    operation_cost::Float64
end

mutable struct Chromosome
    genes::Vector{Int64}
    weighted_time::Float64
    fitness::Float64
    power::Float64
    total_cost::Float64
    depots::Vector{Int64}
    tours::Vector{Tour}
    reds_path::Vector{Tour}
    r_cost::Float64
    r_obj::Float64
    # feasible::Bool # We will consider the infeasible chromosomes later
end


function gene_to_route(gene::Vector{Int64})

    gene_temp = copy(gene)
    deleteat!(gene_temp, findall(x -> x == -1, gene_temp))

    routes = Dict()
    flag = true
    depot = -1
    for i in eachindex(gene_temp)
        if flag
            depot = gene_temp[i]
            if haskey(routes, depot) == false
                routes[depot] = Vector{Int64}()
            end
            flag = false
            continue
        end
        if gene_temp[i] == depot
            flag = true
            continue
        end
        push!(routes[depot], gene_temp[i])
    end

    return routes, collect(keys(routes))
end


# function parent_selection_RWS(population::Vector{Chromosome}, total_p::Float64, popsize::Int64)  #Random Wheel Selection
#     r = -rand() * total_p
#     summ = 0
#     for i in 1:popsize
#         summ -= population[i].power
#         if r > summ
#             return population[i]
#         end
#     end
# end

function parent_selection_TS(population::Vector{Chromosome}, k::Int64, popsize::Int64)  #Tournament Selection
    idx = sample(1:popsize, k, replace=false)
    return population[idx[argmin(idx)]]
end


function parent_selection_RkS(population::Vector{Chromosome}, ranks::Vector{Int64}, total_r::Int64, popsize::Int64)  #Rank Selection
    r = rand() * total_r
    summ = 0
    for i in 1:popsize
        summ += ranks[i]
        if r < summ
            return population[i]
        end
    end
end

function select_parents(population::Vector{Chromosome}, k_tournament::Int64, popsize::Int64)
    return parent_selection_TS(population, k_tournament, popsize), parent_selection_TS(population, k_tournament, popsize)
end

function reproduce(TT::Matrix{Float64}, parent1::Chromosome, parent2::Chromosome, crossover_functions::Vector{Int}, F::Vector{Int64}, G::Vector{Int64}, R::Vector{Int64}, H::Vector{Int64}, 
    distance::Matrix{Float64}, service::Vector{Int64}, patient::Vector{Int64}, demand::Vector{Int64}, capacity_v::Int64, capacity_room::Vector{Int64}, capacity_man::Vector{Int64}, FC::Int64, EC_G::Int64, VC_G::Int64, EC_R::Int64, VC_R::Int64,
    WG::Int64, WR::Int64, epsilon::Int64)

    n_nodes = length(G) + length(F)
    n_depot = length(F)
    
    child = tour_crossover2(parent1, parent2, TT, n_nodes, n_depot, distance, service, patient, demand, capacity_v, EC_G, VC_G)
    using_depots = check_depot(child)
    r_paths, r_cost, r_obj = exact_assignment_algorithm_for_red_(using_depots, R, H, TT, distance, service, patient, demand, capacity_room, capacity_man, VC_R, EC_R, n_nodes)

    obj = child[1].time
    chrm_cost = calculate_costs(child, r_paths, distance, FC, EC_G, VC_G) + r_cost
    S = change_tour_to_gene(child)

    if chrm_cost > epsilon
        if parent1.weighted_time <= parent2.weighted_time
            return parent1
        end
        if parent1.weighted_time > parent2.weighted_time 
            return parent2 
        end
    end

    return Chromosome(S, WG*obj+WR*r_obj, obj, 0.0, chrm_cost, collect(union(Set(using_depots), Set(check_depot(r_paths)))), child, r_paths, r_cost, r_obj)

end

function find_difference(c1::Vector{Int64}, c2::Vector{Int64}, chrm1::Chromosome, chrm2::Chromosome, F::Vector{Int64})  #range between zero and 1, zero when two chromosomes are exactly the same, 1 when all genes are different
    cc1, cc2 = copy(c1), copy(c2)

    for i in vcat(-1, F) deleteat!(cc1, findall(x -> x == i, cc1)) end
    for i in vcat(-1, F) deleteat!(cc2, findall(x -> x == i, cc2)) end

    diff1 = 0
    diff2 = 0

    cc3 = reverse(cc2)
    for i in 1:length(cc1)
        if cc1[i] != cc2[i]
            diff1 += 1
        end
        if cc1[i] != cc3[i]
            diff2 += 1
        end
    end
    return min(diff1, diff2) / length(cc1)
end

function find_difference(c1::Chromosome, c2::Chromosome)  #range between zero and 1, zero when two chromosomes are exactly the same, 1 when all genes are different
    m = length(c1.tours)
    n = length(c1.genes)
    A = zeros(Int, m, m)
    for i in 1:m
        for j in 1:m
            A[i, j] = length(intersect(Set(c1.tours[i].sequence), Set(c2.tours[j].sequence)))
        end
    end
    summ = 0
    while true
        idx = argmax(A)

        if A[idx] == 0
            break
        end
        i = idx[1]
        j = idx[2]
        summ += A[i, j]
        A[i, :] .= 0
        A[:, j] .= 0
    end
    return 1 - summ / n
end

function find_neighbors(a::Int, i::Int, m::Int)
    neighbors = Int[]
    for j = max(1, i - m):i-1
        push!(neighbors, j)
    end
    for j = i+1:min(a, i + m)
        push!(neighbors, j)
    end
    return neighbors
end

function sort_based_on_power!(population::Vector{Chromosome}, num_nei::Int, F::Vector{Int64})
    popsize = length(population)
    diff1 = 0.0
    diff2 = 0.0
    for i in 1:popsize
        neighbors = find_neighbors(popsize, i, num_nei)
        diff = 0.0
        for j in neighbors
            diff1 = find_difference(population[i].genes, population[j].genes, population[i], population[j], F)
            diff += diff1
        end
        population[i].power = population[i].weighted_time * 0.8^(diff / length(neighbors))

    end
    sort!(population, by=x -> x.power)
end


function perform_survival_plan!(population::Vector{Chromosome}, mu::Int64, sigma::Int64)
    if length(population) >= mu + sigma
        del_count = 0
        del_idx = Vector{Int64}()
        @inbounds for i in 1:length(population)-1
            if del_count == length(population) - mu
                break
            end
            @inbounds for j = i+1:length(population)
                if del_count == length(population) - mu
                    break
                end
                if population[i].genes == population[j].genes
                    if !(j in del_idx)
                        push!(del_idx, j)
                        del_count += 1
                    end
                end
            end
        end
        deleteat!(population, sort(del_idx))
        del_idx = Vector{Int64}()
        last_index = length(population)
        index = 0

        @inbounds while del_count < sigma
            i = last_index - index
            push!(del_idx, i)
            del_count += 1
            index += 1
        end
        deleteat!(population, sort(del_idx))
    end
end


function best_route(chrm_::Chromosome)
    for tour in chrm_.tours
        for i in tour.sequence
            print(i, ", ")
        end
        println()
    end
end

function educate_and_add_the_offspring!(offspring::Chromosome, population::Vector{Chromosome}, TT::Matrix{Float64}, Close_nodes::Matrix{Bool}, Customers::Matrix{Float64}, 
    Depots::Matrix{Float64}, old_best::Float64, roullet::Vector{Int}, n_nodes::Int, improve_count::Int, distance::Matrix{Float64}, service::Vector{Int64}, patient::Vector{Int64}, 
    demand::Vector{Int64}, capacity_v::Int64, FC::Int64, EC_G::Int64, VC_G::Int64, WG::Int64, WR::Int64)
    
    if n_nodes > 700 && length(offspring.tours) < 10
        if rand() < 0.1 && improve_count > 100
            solve_all_intersections!(offspring, Customers, Depots, TT, service, patient, demand, capacity_v)
        end
    else
        if rand() < 0.3 
            solve_all_intersections!(offspring, Customers, Depots, TT, service, patient, demand, capacity_v)
            enrich_the_chromosome2!(offspring, TT, Customers, Depots, n_nodes, service, patient, demand, capacity_v)
        end
        if rand() < 0.3 
            solve_all_intersections!(offspring, Customers, Depots, TT, service, patient, demand, capacity_v)
            enrich_the_chromosome2!(offspring, TT, Customers, Depots, n_nodes, service, patient, demand, capacity_v)
        end
    end

    # println(offspring)

    enrich_the_chromosome2!(offspring, TT, Customers, Depots, n_nodes, service, patient, demand, capacity_v)

    # println(offspring)
    if improve_count > 100
        Improve_chromosome!(offspring, TT, Close_nodes, n_nodes, roullet, old_best, 1000, service, patient, demand, capacity_v)
    else
        Improve_chromosome!(offspring, TT, Close_nodes, n_nodes, roullet, old_best, 100, service, patient, demand, capacity_v)
    end
    # println("done educate_and_add_the_offspring!")

    new_genes = Int[]
    delete_tours = Int[]
    obj = 0.0
    for (idx, tour) in enumerate(offspring.tours)
        if tour.staff == 0
            push!(delete_tours, idx)
            continue
        end
        
        for city in tour.sequence
            push!(new_genes, city)
        end
        push!(new_genes, tour.sequence[1])
        push!(new_genes, -1)

        tour.time = find_tour_length(tour.sequence[2:end], TT, tour.sequence[1], service)

        if obj < tour.time
            obj = tour.time
        end

        tour.operation_cost = 0.0
    end
    deleteat!(offspring.tours, delete_tours)
    offspring.genes = new_genes
    offspring.weighted_time = WG*obj+WR*offspring.r_obj
    offspring.fitness = obj
    offspring.total_cost = calculate_costs(offspring.tours, offspring.reds_path, distance, FC, EC_G, VC_G) + offspring.r_cost

    # println(offspring)

    push!(population, offspring)
end

function generate_new_generation(node_size::Int64, TT::Matrix{Float64}, close_nodes::Matrix{Bool}, population::Vector{Chromosome}, popsize::Tuple{Int64,Int64}, 
    k_tournament::Int64, gen_num::Int64, old_best::Float64, old_best_P::Chromosome, improve_count::Int64, mutation_chance::Float64, roullet::Vector{Int}, num_nei::Int, crossover_functions::Vector{Int},
    customers_::Matrix{Float64}, depots_::Matrix{Float64}, F::Vector{Int64}, G::Vector{Int64}, R::Vector{Int64}, H::Vector{Int64}, distance::Matrix{Float64}, service::Vector{Int64}, 
    patient::Vector{Int64}, demand::Vector{Int64}, capacity_vg::Int64, capacity_vr::Int64, capacity_room::Vector{Int64}, capacity_man::Vector{Int64}, FC::Int64, EC_G::Int64, EC_R::Int64, VC_G::Int64, VC_R::Int64, 
    WG::Int64, WR::Int64, coordinate::Vector{Vector}, t0::Float64, time_limit::Float64, epsilon::Int64, obj_vec::Vector{Float64}, verbose::Bool)

    t1 = time()

    mu, sigma = popsize
    n_nodes = node_size 

    if improve_count % 50 == 49
        population = rebalancing(population, TT, customers_, depots_, F, G, R, H, distance, service, patient, demand, capacity_vg, capacity_vr, capacity_room, capacity_man, 
                             FC, EC_G, EC_R, VC_G, VC_R, WG, WR, coordinate, epsilon)
    end

    # if improve_count % 1000 == 999
    #     diversify!(population, TT, K, mu, tsp_tours, customers_, depot_vec, opening_depot, coordinate, F, G, improve_count)
    # end

    # sort_based_on_power!(population, num_nei, F)

    sort!(population, by=x -> x.weighted_time)

    psize = length(population)

    if rand() < mutation_chance
        offspring = mutate(population[rand(1:5)], TT, distance, service, patient, demand, capacity_vg, FC, EC_G, VC_G, WG, WR, epsilon)
    else
        parent1, parent2 = select_parents(population, k_tournament, psize) 
        offspring = reproduce(TT, parent1, parent2, crossover_functions, F, G, R, H, distance, service, patient, demand, capacity_vg, capacity_room, capacity_man, 
                              FC, EC_G, VC_G, EC_R, VC_R, WG, WR, epsilon)
    end

    educate_and_add_the_offspring!(offspring, population, TT, close_nodes, customers_, depots_, old_best, roullet, n_nodes, improve_count, 
                                   distance, service, patient, demand, capacity_vg, FC, EC_G, VC_G, WG, WR)

    sort!(population, by=x -> x.weighted_time)

    perform_survival_plan!(population, mu, sigma)

    new_best = population[1].weighted_time
    new_best_g_obj = population[1].fitness
    new_best_r_obj = population[1].r_obj
    new_best_P = population[1]

    if round(old_best, digits=3) > round(new_best, digits=3)
        push!(obj_vec, new_best)
        old_best = new_best
        old_best_P = new_best_P
        improve_count = 0
    else
        push!(obj_vec, old_best)
        improve_count += 1
    end
    t2 = time()

    if verbose
        if gen_num % 100 == 0
            println("Generation ", gen_num, " the best objective is: ", old_best, "   time left: $(round(t0+time_limit -time())) seconds")
        end
    end
    gen_num += 1

    return gen_num, old_best, old_best_P, improve_count, obj_vec
end

function perform_genetic_algorithm(
    node_size::Int64, TT::Matrix{Float64}, h::Float64, popsize::Tuple{Int64,Int64}, k_tournament::Int64, num_iter::Int64, time_limit::Float64, mutation_chance::Float64, 
    num_nei::Int64, crossover_functions::Vector{Int}, Customers_::Matrix{Float64}, Depots_::Matrix{Float64}, F::Vector{Int64}, G::Vector{Int64}, R::Vector{Int64}, H::Vector{Int64}, 
    Return_Depot::Vector{Int64}, distance::Matrix{Float64}, service::Vector{Int64}, patient::Vector{Int64}, demand::Vector{Int64}, capacity_vg::Int64, capacity_vr::Int64, capacity_room::Vector{Int64}, capacity_man::Vector{Int64}, 
    FC::Int64, EC_G::Int64, EC_R::Int64, VC_G::Int64, VC_R::Int64, WG::Int64, WR::Int64, coordinate::Vector{Vector}, epsilon::Int64, verbose::Bool)    

    Random.seed!(Int(round(time())))
    t1 = time()
    
    ClosenessT = find_closeness(TT, length(F), length(G), h) # return number of (customer node)*h closeness node to node time
    mu, sigma = popsize
    improve_count = 0
    Gen_num = 0
    old_best = 0.0
    roullet = ones(Int, 4) * 100

    depot_vec, num_depot_vec = open_depot(F, mu)
    assignments = assign_patient(depot_vec, Customers_)
    tours = allocate_depot(F, assignments, depot_vec, coordinate)
    tsp_tours, using_depot = find_mdtsp_tour(tours, distance, coordinate)
    r_obj, r_cost, r_tour = exact_assignment_algorithm_for_red(using_depot, F, G, R, H, TT, distance, service, patient, demand, capacity_room, capacity_man, VC_R, EC_R)
    Population, old_best = generate_initial_population(TT, mu, sigma, num_depot_vec, tsp_tours, r_obj, r_cost, r_tour, Customers_, F, G, R, H, Return_Depot, 
                                                       distance, service, patient, demand, capacity_vg, capacity_vr, FC, EC_G, EC_R, VC_G, VC_R,
                                                       WG, WR)
    old_best_P = Population[1]

    if verbose
        println("The initialization took ", time() - t1, " seconds.")
        println("The initial obj is ", old_best)
        println("The initial G obj is ", Population[1].fitness)
        println("The initial R obj is ", Population[1].r_obj)
        println("The initial Operation cost is ", Population[1].total_cost)
    end

    obj_vec = Vector{Float64}()
    push!(obj_vec, old_best)
    while improve_count < num_iter
        if time() - t1 >= time_limit
            break
        end

        Gen_num, old_best, old_best_P, improve_count, obj_vec = generate_new_generation(node_size, TT, ClosenessT, Population, popsize, k_tournament, Gen_num, old_best, old_best_P,
                                                                               improve_count, mutation_chance, roullet, num_nei, crossover_functions, Customers_, 
                                                                               Depots_, F, G, R, H, distance, service, patient, demand, capacity_vg, capacity_vr, capacity_room, capacity_man, 
                                                                               FC, EC_G, EC_R, VC_G, VC_R, WG, WR, coordinate, t1, time_limit, epsilon, obj_vec, verbose)
    end
    t2 = time()

    if verbose 
        println("The best objective achieved in ", Gen_num, " generations is: ", old_best_P.weighted_time, " and it took ", t2 - t1, " seconds.")
        println("The G obj is ", old_best_P.fitness)
        println("The R obj is ", old_best_P.r_obj)
        println("The operation cost is ", old_best_P.total_cost, '$')
        println("And the best route is: ")

        best_route(old_best_P)
    end

    return old_best_P, roullet, obj_vec
end