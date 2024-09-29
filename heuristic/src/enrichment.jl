function enrich_the_chromosome2!(chrm::Chromosome, T::Matrix{Float64}, customers::Matrix{Float64}, depots::Matrix{Float64}, n_nodes::Int,
    service::Vector{Int64}, patient::Vector{Int64}, demand::Vector{Int64}, capacity_v::Int64)

    # println("Start enrich_the_chromosome2!")
    m = length(chrm.tours)
    temp = deepcopy(chrm)
    max_tour_index = argmax([chrm.tours[i].time for i in 1:length(chrm.tours)])
    max_tour_length = chrm.fitness
    tour_neighbors = find_tour_neighbors(chrm.tours, T, m)
    improved = true
    count = 0

    while improved && count < 100
        count += 1
        improved = false
        not_swapped = fill(true, m, m)
        for r1 in 1:m
            for r2 in tour_neighbors[r1]
                if r2 != max_tour_index # if max time tour than shift
                    improved = shift_one!(T, chrm, r1, r2, max_tour_index, max_tour_length, improved, service, patient, demand, capacity_v)
                    max_tour_index = argmax([chrm.tours[i].time for i in 1:length(chrm.tours)])
                    max_tour_length = maximum([chrm.tours[i].time for i in 1:length(chrm.tours)])
                end
                if not_swapped[r1, r2] # else swap
                    improved = swap_one!(T, chrm.tours, r1, r2, max_tour_index, max_tour_length, improved, service, patient, demand, capacity_v) 
                    max_tour_index = argmax([chrm.tours[i].time for i in 1:length(chrm.tours)])
                    max_tour_length = maximum([chrm.tours[i].time for i in 1:length(chrm.tours)])
                    not_swapped[r1, r2] = false
                    not_swapped[r2, r1] = false
                end
            end
        end

    end
    for tour in chrm.tours
        two_opt_on_route(tour, T, n_nodes)
    end
    chrm.genes = Int[]
    chrm.fitness = maximum([chrm.tours[i].time for i in 1:length(chrm.tours)])
    for tour in chrm.tours
        chrm.genes = vcat(chrm.genes, tour.sequence, tour.sequence[1], -1)
        # tour.operation_cost = 0.0
    end
    # chrm.total_cost = calculate_costs(chrm.tours, distance, FC, EC_G, VC_G) + r_cost
end


function shift_one!(T::Matrix{Float64}, chrm::Chromosome, r1::Int, r2::Int, max_tour_index::Int, max_tour_length::Float64, improved::Bool,
    service::Vector{Int64}, patient::Vector{Int64}, demand::Vector{Int64}, capacity_v::Int64)

    if chrm.tours[r1].time <= chrm.tours[r2].time
        return improved
    end
    tour1 = chrm.tours[r1].sequence[2:end]
    depot1 = chrm.tours[r1].sequence[1]
    tour2 = chrm.tours[r2].sequence[2:end]
    depot2 = chrm.tours[r2].sequence[1]

    k1 = 0
    while k1 < length(tour1)
        time1 = chrm.tours[r1].time
        time2 = chrm.tours[r2].time
        k1 += 1
        city1 = tour1[k1]

        if capacity_v - max(chrm.tours[r2].staff, demand[city1]) < sum(chrm.tours[r2].patient) + patient[city1] continue end

        for k2 in 1:length(tour2)+1
            new_time2 = calculate_new_cost_add_one(tour2, depot2, time2, city1, k2, T) + service[city1]
            new_time1 = calculate_new_cost_remove_one(tour1, depot1, time1, k1, T) - service[city1]
            do_it = false
            if r1 == max_tour_index
                if new_time2 < max_tour_length
                    do_it = true
                end
            else
                if (new_time2 - time2) < (time1 - new_time1) && new_time2 < max_tour_length
                    do_it = true
                end
            end
            if do_it
                insert!(tour2, k2, city1)
                deleteat!(tour1, k1)

                chrm.tours[r1].time = new_time1
                chrm.tours[r1].sequence = vcat(depot1, tour1)
                chrm.tours[r1].staff = find_tour_staff(chrm.tours[r1].sequence, demand)
                chrm.tours[r1].patient = find_tour_patient(chrm.tours[r1].sequence[2:end], patient)
                
                chrm.tours[r2].time = new_time2
                chrm.tours[r2].sequence = vcat(depot2, tour2)
                chrm.tours[r2].staff = find_tour_staff(chrm.tours[r2].sequence, demand)
                chrm.tours[r2].patient = find_tour_patient(chrm.tours[r2].sequence[2:end], patient)

                k1 -= 1
                improved = true
                break
            end
        end
    end
    return improved
end

function swap_one!(T::Matrix{Float64}, tours::Vector{Tour}, r1::Int, r2::Int, max_tour_index::Int, max_tour_length::Float64, improved::Bool,
    service::Vector{Int64}, patient::Vector{Int64}, demand::Vector{Int64}, capacity_v::Int64)

    tour1 = tours[r1].sequence[2:end]
    depot1 = tours[r1].sequence[1]
    tour2 = tours[r2].sequence[2:end]
    depot2 = tours[r2].sequence[1]
    k1 = 0
    while k1 < length(tour1)
        time1 = tours[r1].time
        time2 = tours[r2].time
        k1 += 1
        city1 = tour1[k1]
        for k2 in 1:length(tour2)
            city2 = tour2[k2]
            temp1 = replace(tour1, city1=>city2)
            temp2 = replace(tour2, city2=>city1)

            if capacity_v - find_tour_staff(temp2, demand) < sum(temp2) continue end
            if capacity_v - find_tour_staff(temp1, demand) < sum(temp1) continue end

            new_time1, new_time2 = calculate_new_cost_swap_one(tour1, depot1, time1, city1, k1, tour2, depot2, time2, city2, k2, T)
            new_time1 = new_time1 - service[city1] + service[city2]
            new_time2 = new_time2 - service[city2] + service[city1]

            do_it = false
            if r1 == max_tour_index
                if new_time1 < time1 && new_time2 < time1
                    do_it = true
                end
            elseif r2 == max_tour_index
                if new_time1 < time2 && new_time2 < time2
                    do_it = true
                end
            else
#                 if (new_cost1 + new_cost2 < cost1 + cost2) && new_cost2 < max_tour_length && new_cost1 < max_tour_length
                if new_time2 < time2 && new_time1 < time1
                    do_it = true
                end
            end
            if do_it
                tours[r1].time = new_time1
                tours[r2].time = new_time2

                tour1[k1] = city2
                tour2[k2] = city1

                tours[r1].sequence = vcat(depot1, tour1)
                tours[r1].staff = find_tour_staff(tours[r1].sequence, demand)
                tours[r1].patient = find_tour_patient(tours[r1].sequence[2:end], patient)
                
                tours[r2].sequence = vcat(depot2, tour2)
                tours[r2].staff = find_tour_staff(tours[r2].sequence, demand)
                tours[r2].patient = find_tour_patient(tours[r2].sequence[2:end], patient)
                
                improved = true
                break
            end
        end
    end
    return improved
end

function shift_two!(T::Matrix{Float64}, tours::Vector{Tour}, r1::Int, r2::Int, max_tour_index::Int, max_tour_length::Float64, n_nodes::Int, improved::Bool) 
    if tours[r1].cost <= tours[r2].cost
        return improved
    end
    tour1 = tours[r1].sequence
    tour2 = tours[r2].sequence
    k1 = 0
    while k1 < length(tour1)-1
        cost1 = tours[r1].cost
        cost2 = tours[r2].cost
        k1 += 1
        city1 = tour1[k1]
        city2 = tour1[k1+1]
        for k2 in 1:length(tour2)+1
            
            new_cost2, straight = calculate_new_cost_add_two(tour2, cost2, city1, city2, k2, T, n_nodes)
            new_cost1 = calculate_new_cost_remove_two(tour1, cost1, k1, T, n_nodes)
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
                tours[r1].cost = new_cost1
                tours[r2].cost = new_cost2
                if straight
                    insert!(tour2, k2, city1)
                    insert!(tour2, k2 + 1, city2)
                else
                    insert!(tour2, k2, city2)
                    insert!(tour2, k2 + 1, city1)
                end
                deleteat!(tour1, [k1, k1 + 1])
                k1 -= 1
                improved = true
                break
            end
        end
    end
    return improved
end

function swap_two!(T::Matrix{Float64}, tours::Vector{Tour}, r1::Int, r2::Int, max_tour_index::Int, max_tour_length::Float64, n_nodes::Int, improved::Bool)    
    tour1 = tours[r1].sequence
    tour2 = tours[r2].sequence
    k1 = 0
    while k1 < length(tour1)-1
        cost1 = tours[r1].cost
        cost2 = tours[r2].cost
        k1 += 1
        city11 = tour1[k1]
        city12 = tour1[k1+1]
        for k2  in 1:length(tour2)-1
            city21 = tour2[k2]
            city22 = tour2[k2+1]
            new_cost1, new_cost2, straight1, straight2 = calculate_new_cost_swap_two_updated(tour1, cost1, city11, city12, k1, tour2, cost2, city21, city22, k2, T, n_nodes)
            do_it = false
            if r1 == max_tour_index
                if new_cost1 < cost1 && new_cost2 < cost1
                    do_it = true
                end
            elseif r2 == max_tour_index
                if new_cost1 < cost2 && new_cost2 < cost2
                    do_it = true
                end
            else
#                 if (new_cost1 + new_cost2 < cost1 + cost2) && new_cost2 < max_tour_length && new_cost1 < max_tour_length
                if new_cost2 < cost2 && new_cost1 < cost1
                    do_it = true
                end
            end
            if do_it
                tours[r1].cost = new_cost1
                tours[r2].cost = new_cost2
                if straight2
                    tour1[k1] = city21
                    tour1[k1+1] = city22
                else
                    tour1[k1] = city22
                    tour1[k1+1] = city21
                end
                if straight1
                    tour2[k2] = city11
                    tour2[k2+1] = city12
                else
                    tour2[k2] = city12
                    tour2[k2+1] = city11
                end
                improved = true
                break
            end
        end
    end
    return improved
end

function cross_two_points!(T::Matrix{Float64}, tours::Vector{Tour}, r1::Int, r2::Int, max_tour_index::Int, max_tour_length::Float64, n_nodes::Int, improved::Bool)    
    tour1 = tours[r1].sequence
    tour2 = tours[r2].sequence
    if length(tour1) < 6 || length(tour2) < 6
        return improved
    end
    k11 = 0
    while k11 < length(tour1)-5
        cost1 = tours[r1].cost
        cost2 = tours[r2].cost
        k11 += 1
        for tau1 in 3:5
            k12 = k11 + tau1
            for k21 in 1:length(tour2)-5
                for tau2 in 3:5
                    k22 = k21 + tau2
                    new_cost1, new_cost2, straight1, straight2 = calculate_new_cost_cross(tour1, cost1, tour2, cost2, k11, k12, k21, k22, T, n_nodes)
                    do_it = false
                    if r1 == max_tour_index
                        if new_cost1 < cost1 && new_cost2 < cost1
                            do_it = true
                        end
                    elseif r2 == max_tour_index
                        if new_cost1 < cost2 && new_cost2 < cost2
                            do_it = true
                        end
                    else
        #                 if (new_cost1 + new_cost2 < cost1 + cost2) && new_cost2 < max_tour_length && new_cost1 < max_tour_length
                        if new_cost2 < cost2 && new_cost1 < cost1
                            do_it = true
                        end
                    end
                    if do_it
                        tours[r1].cost = new_cost1
                        tours[r2].cost = new_cost2
                        if straight2
                            alpha1 = copy(tour1[k11:k12])
                        else
                            alpha1 = reverse(copy(tour1[k11:k12]))
                        end
                        if straight1
                            alpha2 = copy(tour2[k21:k22])
                        else
                            alpha2 = reverse(copy(tour2[k21:k22]))
                        end
                        deleteat!(tour1, [i for i = k11:k12])
                        for i in 1:k22-k21+1
                            insert!(tour1, i + k11 - 1, alpha2[i])
                        end

                        deleteat!(tour2, [i for i = k21:k22])
                        for i in 1:k12-k11+1
                            insert!(tour2, i + k21 - 1, alpha1[i])
                        end
#                         println("Done")
                        return true
                    end
                end
            end
        end
    end
    return improved
end