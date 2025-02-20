function crossover_POS(parent1::Vector{Int64}, parent2::Vector{Int64}, n_nodes::Int64)  #position based crossover
    child = zeros(Int64, n_nodes)
    num_pos = rand(1:n_nodes-1)
    selected_pos = sample(1:n_nodes, Weights(ones(n_nodes)), num_pos, replace=false)
    selected_p1 = parent1[selected_pos]
    child[selected_pos] = selected_p1

    for i in parent2
        if !(i in selected_p1)
            child[findfirst(x -> x == 0, child)] = i
        end
    end
    return child
end


function crossover_OX2(parent1::Vector{Int64}, parent2::Vector{Int64}, n_nodes::Int64)   #Order-based crossover
    child = zeros(Int64, n_nodes)
    num_pos = rand(1:n_nodes)
    selected_pos2 = sample(1:n_nodes, Weights(ones(n_nodes)), num_pos, replace=false)
    selected_p2 = parent2[selected_pos2]

    selected_pos1 = findall(x -> x in selected_p2, parent1)
    unselected_pos = setdiff(1:n_nodes, selected_pos1)
    child[unselected_pos] = parent1[unselected_pos]
    child[selected_pos1] = parent2[sort(selected_pos2)]
    return child
end

function crossover_HX_(TT::Matrix{Float64}, parent1::Vector{Int64}, parent2::Vector{Int64}, n_nodes::Int64)   #heuristic crossover 

    remaining_cities = copy(parent1)
    r = rand(1:n_nodes)
    current_city = remaining_cities[r]
    child = Int[current_city]
    deleteat!(remaining_cities, r)
    p = 1
    child_mask = [false for i in 1:n_nodes]
    child_mask[current_city] = true

    while sum(child_mask) < n_nodes
        if p == 1
            pos1 = r
            pos2 = findfirst(x -> x == current_city, parent2)
        else
            pos1 = findfirst(x -> x == current_city, parent1)
            pos2 = r
        end

        next_city = n_nodes + 2
        min_edge = 100000

        if pos1 == 1
            if !child_mask[parent1[pos1+1]]
                if TT[current_city+1, parent1[pos1+1]+1] < min_edge
                    next_city = parent1[pos1+1]
                    min_edge = TT[current_city+1, parent1[pos1+1]+1]
                    p = 1
                    r = pos1 + 1
                end
            end
        elseif pos1 == n_nodes
            if !child_mask[parent1[pos1-1]]
                if TT[current_city+1, parent1[pos1-1]+1] < min_edge
                    next_city = parent1[pos1-1]
                    min_edge = TT[current_city+1, parent1[pos1-1]+1]
                    p = 1
                    r = pos1 - 1
                end
            end
        else
            if !child_mask[parent1[pos1+1]]
                if TT[current_city+1, parent1[pos1+1]+1] < min_edge
                    next_city = parent1[pos1+1]
                    min_edge = TT[current_city+1, parent1[pos1+1]+1]
                    p = 1
                    r = pos1 + 1
                end
            end
            if !child_mask[parent1[pos1-1]]
                if TT[current_city+1, parent1[pos1-1]+1] < min_edge
                    next_city = parent1[pos1-1]
                    min_edge = TT[current_city+1, parent1[pos1-1]+1]
                    p = 1
                    r = pos1 - 1
                end
            end
        end

        if pos2 == 1
            if !child_mask[parent2[pos2+1]]
                if TT[current_city+1, parent2[pos2+1]+1] < min_edge
                    next_city = parent2[pos2+1]
                    min_edge = TT[current_city+1, parent2[pos2+1]+1]
                    p = 2
                    r = pos2 + 1
                end
            end
        elseif pos2 == n_nodes
            if !child_mask[parent2[pos2-1]]
                if TT[current_city+1, parent2[pos2-1]+1] < min_edge
                    next_city = parent2[pos2-1]
                    min_edge = TT[current_city+1, parent2[pos2-1]+1]
                    p = 2
                    r = pos2 - 1
                end
            end
        else
            if !child_mask[parent2[pos2+1]]
                if TT[current_city+1, parent2[pos2+1]+1] < min_edge
                    next_city = parent2[pos2+1]
                    min_edge = TT[current_city+1, parent2[pos2+1]+1]
                    p = 2
                    r = pos2 + 1
                end
            end
            if !child_mask[parent2[pos2-1]]
                if TT[current_city+1, parent2[pos2-1]+1] < min_edge
                    next_city = parent2[pos2-1]
                    min_edge = TT[current_city+1, parent2[pos2-1]+1]
                    p = 2
                    r = pos2 - 1
                end
            end
        end

        if next_city == n_nodes + 2
            remainings = findall(x -> x == false, child_mask)
            next_city = remainings[rand(1:length(remainings))]
        end

        current_city = next_city
        push!(child, next_city)
        child_mask[next_city] = true
    end
    return child
end

function crossover_HX(TT::Matrix{Float64}, parent1::Vector{Int64}, parent2::Vector{Int64}, n_nodes::Int64)   #heuristic crossover 
    # 기존 input 되는 parent는 customer node로만 구성

    remaining_cities = copy(parent1)
    current_city = rand(parent1)
    child = Int[current_city]
    deleteat!(remaining_cities, findfirst(x -> x == current_city, remaining_cities))
    while length(remaining_cities) > 0
        pos1 = findfirst(x -> x == current_city, parent1)
        pos2 = findfirst(x -> x == current_city, parent2)
        next_city = n_nodes + 2
        min_edge = 100000

        if pos1 == 1 # change to first node
            if !(parent1[pos1] in child)
                if TT[current_city, parent1[pos1+1]] < min_edge
                    next_city = parent1[pos1+1]
                    min_edge = TT[current_city, parent1[pos1+1]]
                end
            end
        elseif pos1 == n_nodes
            if !(parent1[pos1-1] in child)
                if TT[current_city, parent1[pos1-1]] < min_edge
                    next_city = parent1[pos1-1]
                    min_edge = TT[current_city, parent1[pos1-1]]
                end
            end
        else
            if !(parent1[pos1+1] in child)
                if TT[current_city, parent1[pos1+1]] < min_edge
                    next_city = parent1[pos1+1]
                    min_edge = TT[current_city, parent1[pos1+1]]
                end
            end
            if !(parent1[pos1-1] in child)
                if TT[current_city, parent1[pos1-1]] < min_edge
                    next_city = parent1[pos1-1]
                    min_edge = TT[current_city, parent1[pos1-1]]
                end
            end
        end

        if pos2 == 1
            if !(parent2[pos2+1] in child)
                if TT[current_city, parent2[pos2+1]] < min_edge
                    next_city = parent2[pos2+1]
                    min_edge = TT[current_city, parent2[pos2+1]]
                end
            end
        elseif pos2 == n_nodes
            if !(parent2[pos2-1] in child)
                if TT[current_city, parent2[pos2-1]] < min_edge
                    next_city = parent2[pos2-1]
                    min_edge = TT[current_city, parent2[pos2-1]]
                end
            end
        else
            if !(parent2[pos2+1] in child)
                if TT[current_city, parent2[pos2+1]] < min_edge
                    next_city = parent2[pos2+1]
                    min_edge = TT[current_city, parent2[pos2+1]]
                end
            end
            if !(parent2[pos2-1] in child)
                if TT[current_city, parent2[pos2-1]] < min_edge
                    next_city = parent2[pos2-1]
                    min_edge = TT[current_city, parent2[pos2-1]]
                end
            end
        end

        # If do not find (current to next) shorter than min_edge, randomly select next city
        if next_city == n_nodes + 2
            next_city = remaining_cities[rand(1:length(remaining_cities))]
        end

        current_city = next_city
        push!(child, next_city)
        deleteat!(remaining_cities, findfirst(x -> x == current_city, remaining_cities))
    end
    return child
end

function remove_cities_from_one_tour(tour_::Tour, cities::Vector{Int}, T::Matrix{Float64})
    # println("Start remove_cities_from_one_tour")
    nt = length(tour_.sequence[2:end])
    index = 1
    i = 1
    seq = Int[]
    seqs = Vector{Vector{Int}}()
    while i <= cities[length(cities)]
        if i == cities[index]
            push!(seq, i)
            if i == cities[length(cities)]
                push!(seqs, seq)
            end
            i += 1
            index += 1
        else
            if length(seq) > 0
                push!(seqs, seq)
                seq = Int[]
            end
            i += 1
        end
    end
    tour = tour_.sequence[2:end]
    depot = tour_.sequence[1]
    time = tour_.time
    
    for seq in seqs
        ns = length(seq)
        if seq[1] == 1
            if ns == 1
                time = time - T[depot, tour[1]] - T[tour[1], tour[2]] + T[depot, tour[2]]
            else
                time = time - T[depot, tour[1]] - sum(T[tour[seq[i]], tour[seq[i+1]]] for i in 1:ns-1) - T[tour[seq[ns]], tour[seq[ns]+1]] + T[depot, tour[seq[ns]+1]]
            end
        elseif seq[ns] == nt
            if ns == 1
                time = time - T[tour[nt], depot] - T[tour[nt-1], tour[nt]] + T[tour[nt-1], depot]
            else
                time = time - T[tour[nt], depot] - sum(T[tour[seq[i]], tour[seq[i+1]]] for i in 1:ns-1) - T[tour[seq[1]-1], tour[seq[1]]] + T[tour[seq[1]-1], depot]
            end
        else
            if ns == 1
                time = time - T[tour[seq[1]-1], tour[seq[1]]] - T[tour[seq[1]], tour[seq[1]+1]] + T[tour[seq[1]-1], tour[seq[1]+1]]
            else
                time = time - T[tour[seq[1]-1], tour[seq[1]]] - sum(T[tour[seq[i]], tour[seq[i+1]]] for i in 1:ns-1) - T[tour[seq[ns]], tour[seq[ns]+1]] + T[tour[seq[1]-1], tour[seq[ns]+1]]
            end
        end
    end
    tour_.time = time
end


# function remove_cities_from_one_tour(tour_::Tour, cities::Vector{Int}, T::Matrix{Float64}, n_nodes::Int)
#     # println("Start remove_cities_from_one_tour")
#     nt = length(tour_.sequence[2:end])
#     index = 1
#     i = 1
#     seq = Int[]
#     seqs = Vector{Vector{Int}}()
#     while i <= cities[length(cities)]
#         if i == cities[index]
#             push!(seq, i)
#             if i == cities[length(cities)]
#                 push!(seqs, seq)
#             end
#             i += 1
#             index += 1
#         else
#             if length(seq) > 0
#                 push!(seqs, seq)
#                 seq = Int[]
#             end
#             i += 1
#         end
#     end
#     tour = tour_.sequence[2:end]
#     depot = tour_.sequence[1]
#     time = tour_.time
    
#     for seq in seqs
#         ns = length(seq)
#         if seq[1] == 1
#             if ns == 1
#                 time = time - T[depot, tour[1]] - T[tour[1], tour[2]] + T[depot, tour[2]]
#             else
#                 time = time - T[depot, tour[1]] - sum(T[tour[seq[i]], tour[seq[i+1]]] for i in 1:ns-1) - T[tour[seq[ns]], tour[seq[ns]+1]] + T[depot, tour[seq[ns]+1]]
#             end
#         elseif seq[ns] == nt
#             if ns == 1
#                 time = time - T[tour[nt], depot] - T[tour[nt-1], tour[nt]] + T[tour[nt-1], depot]
#             else
#                 time = time - T[tour[nt], depot] - sum(T[tour[seq[i]], tour[seq[i+1]]] for i in 1:ns-1) - T[tour[seq[1]-1], tour[seq[1]]] + T[tour[seq[1]-1], depot]
#             end
#         else
#             if ns == 1
#                 time = time - T[tour[seq[1]-1], tour[seq[1]]] - T[tour[seq[1]], tour[seq[1]+1]] + T[tour[seq[1]-1], tour[seq[1]+1]]
#             else
#                 time = time - T[tour[seq[1]-1], tour[seq[1]]] - sum(T[tour[seq[i]], tour[seq[i+1]]] for i in 1:ns-1) - T[tour[seq[ns]], tour[seq[ns]+1]] + T[tour[seq[1]-1], tour[seq[ns]+1]]
#             end
#         end
#     end
#     tour_.time = time
# end


function Remove_one_city(tours::Vector{Tour}, city::Int, T::Matrix{Float64}, n_nodes::Int)
    for (t, tour) in enumerate(tours)
        nt = length(tour.sequence)
        k = findfirst(x -> x == city, tour.sequence)
        if !isnothing(k)
            if k == 1
                if nt == 1
                    deleteat!(tours, t)
                else
                    tour.time = tour.time - T[1, city+1] - T[city+1, tour.sequence[2]+1] + T[1, tour.sequence[2]+1]
                    deleteat!(tour.sequence, 1)
                end
            elseif k == nt
                tour.time = tour.time - T[city+1, n_nodes+2] - T[tour.sequence[nt-1]+1, city+1] + T[tour.sequence[nt-1]+1, n_nodes+2]
                deleteat!(tour.sequence, nt)
            else
                tour.time = tour.time - T[tour.sequence[k-1]+1, city+1] - T[city+1, tour.sequence[k+1]+1] + T[tour.sequence[k-1]+1, tour.sequence[k+1]+1]
                deleteat!(tour.sequence, k)
            end
        end
    end
end

function put_city_in_tour(c::Vector{Tour}, city::Int, T::Matrix{Float64}, n_nodes::Int, service::Vector{Int64})
    least_increase = Inf
    if length(c) >= 2
        best_tour = 0
        best_position = 0
        for i = 2:length(c)
            tour = c[i].sequence
            depot = tour[1]
            nt = length(tour)
            if nt == 1
                increase = T[depot, city] + service[city]
                if increase < least_increase
                    least_increase = increase
                    best_tour = i
                    best_position = 2
                end
            else
                increase = T[depot, city] + service[city] + T[city, tour[2]] - T[depot, tour[2]]
                if increase < least_increase
                    least_increase = increase
                    best_tour = i
                    best_position = 2
                end
                for j = 3:nt
                    increase = T[tour[j-1], city] + service[city] + T[city, tour[j]] - T[tour[j-1], tour[j]]
                    if increase < least_increase
                        least_increase = increase
                        best_tour = i
                        best_position = j
                    end
                end
                increase = T[tour[nt], city] + service[city]
                if increase < least_increase
                    least_increase = increase
                    best_tour = i
                    best_position = nt + 1
                end
            end
        end
        insert!(c[best_tour].sequence, best_position, city)
        c[best_tour].time = find_tour_length(c[best_tour].sequence[2:end], T, c[best_tour].sequence[1], service)
        if c[best_tour].time > c[1].time
            temp = deepcopy(c[1])
            c[1] = c[best_tour]
            c[best_tour] = temp
        end
    else
        best_position = 0
        tour = c[1].sequence
        depot = tour[1]
        nt = length(tour)
        if nt == 1
            increase = T[depot, city] + service[city]
            if increase < least_increase
                least_increase = increase
                best_position = 2
            end
        else
            increase = T[depot, city] + service[city] + T[city, tour[2]] - T[depot, tour[2]]
            if increase < least_increase
                least_increase = increase
                best_position = 2
            end
            for j = 3:nt
                increase = T[tour[j-1], city] + service[city] + T[city, tour[j]] - T[tour[j-1], tour[j]]
                if increase < least_increase
                    least_increase = increase
                    best_position = j
                end
            end
            increase = T[tour[nt], city] + service[city]
            if increase < least_increase
                least_increase = increase
                best_position = nt + 1
            end
        end
        insert!(c[1].sequence, best_position, city)
        c[1].time = find_tour_length(c[1].sequence[2:end], T, c[1].sequence[1], service)
    end
end

function put_cities_in_tour(c::Tour, cities::Vector{Int}, T::Matrix{Float64}, n_nodes::Int)
    tour = c.sequence
    depot = tour[1]
    for city in cities
        least_increase = Inf
        best_position = 0
        nt = length(tour)
        if nt == 1
            increase = T[depot, city] + T[city, depot]
            if increase < least_increase
                least_increase = increase
                best_position = 2
            end
        else
            increase = T[depot, city] + T[city, tour[2]] - T[depot, tour[2]]
            if increase < least_increase
                least_increase = increase
                best_position = 2
            end
            for j = 3:nt
                increase = T[tour[j-1], city] + T[city, tour[j]] - T[tour[j-1], tour[j]]
                if increase < least_increase
                    least_increase = increase
                    best_position = j
                end
            end
            increase = T[tour[nt], city] + T[city, depot] - T[tour[nt], depot]
            if increase < least_increase
                least_increase = increase
                best_position = nt + 1
            end
        end
        insert!(tour, best_position, city)
        c.time += least_increase
    end
end

function tour_crossover2(parent1::Chromosome, parent2::Chromosome, T::Matrix{Float64}, n_nodes::Int64, n_depot::Int64, distance::Matrix{Float64},
     service::Vector{Int64}, patient_vec::Vector{Int64}, demand::Vector{Int64}, capacity_v::Int64, EC_G::Int64, VC_G::Int64)
    #2  At each step, select a tour from parent1, and select the tour with maximum mutual cities from parent2. 
    # Conduct a simple two point crossover between them and add it to the new tours.
    # At the end, Delete all the repeating cities from the tours. 
    # All the remaining cities will be placed in the current tours based on a greedy approach (minimum increase)
    # P1 = deepcopy(parent1)
    # P2 = deepcopy(parent2)
    P1_tours = deepcopy(parent1.tours)
    P2_tours = deepcopy(parent2.tours)
    c = Tour[]
    m = min(length(P1_tours), length(P2_tours))

    for i in 1:m
        tour1 = P1_tours[i].sequence
        time1 = P1_tours[i].time
        staff1 = P1_tours[i].staff
        patient1 = P1_tours[i].patient
        operation_cost1 = P1_tours[i].operation_cost
        depot = tour1[1]
        max_intersection = -1
        tour2 = Int[]
        time2 = 0.0
        staff2 = 0
        patient2 = Int[]
        operation_cost2 = 0.0
        r2 = 0
        for j in eachindex(P2_tours)
            inter = length(intersect(tour1, P2_tours[j].sequence))
            if inter > max_intersection
                max_intersection = inter
                tour2 = P2_tours[j].sequence
                time2 = P2_tours[j].time
                staff2 = P2_tours[j].staff
                patient2 = P2_tours[j].patient
                operation_cost2 = P2_tours[j].operation_cost
                r2 = j
            end
        end

        deleteat!(P2_tours, r2) # delete maximum tour of intersection 

        if length(tour1) <= length(tour2)
            if length(tour1) <= 4
                push!(c, Tour(tour2, time2, staff2, patient2, 0.0))
            else
                idx1, idx2 = sort(sample(2:length(tour1)-1, 2, replace=false))
                cc = vcat(depot, tour2[2:idx1-1], tour1[idx1:idx2], tour2[idx2+1:length(tour2)])
                staff = find_tour_staff(cc, demand)
                patient = find_tour_patient(cc[2:end], patient_vec)
                if capacity_v - staff < sum(patient)
                    push!(c, Tour(tour2, time2, staff2, patient2, 0.0))
                else
                    push!(c, Tour(cc, find_tour_length(cc[2:end], T, depot, service), staff, patient, 0.0))
                end
            end
        else
            if length(tour2) <= 4
                push!(c, Tour(tour1, time1, staff1, patient1, 0.0))
            else
                idx1, idx2 = sort(sample(2:length(tour2)-1, 2, replace=false))
                cc = vcat(depot, tour1[2:idx1-1], tour2[idx1:idx2], tour1[idx2+1:length(tour1)])
                staff = find_tour_staff(cc, demand)
                patient = find_tour_patient(cc[2:end], patient_vec)
                if capacity_v - staff < sum(patient)
                    push!(c, Tour(tour1, time1, staff1, patient1, 0.0))
                else
                    push!(c, Tour(cc, find_tour_length(cc[2:end], T, depot, service), staff, patient, 0.0))
                end
            end
        end
    end

    counters = zeros(n_nodes)
    outsiders = Int[]
    for tour in c
        delete_indices = Int[]
        for (j, node) in enumerate(tour.sequence[2:end])
            if counters[node] > 0
                push!(delete_indices, j)
            else
                counters[node] += 1
            end
        end
        if length(delete_indices) == length(tour.sequence[2:end]) # all of tour nodes are duplicated node 
            tour.time = 0.0
            tour.sequence = Int[tour.sequence[1]]
        elseif length(delete_indices) > 0 # subsection duplicated node 
            # tour.time = remove_cities_from_one_tour(tour, delete_indices, T)
            # deleteat!(tour.sequence, delete_indices.+1)
            # tour.time += sum([service[node] for node in tour.sequence[2:end]])
            deleteat!(tour.sequence, delete_indices.+1)
            tour.time = find_tour_length(tour.sequence[2:end], T, tour.sequence[1], service)
            tour.staff = find_tour_staff(tour.sequence, demand)
            tour.patient = find_tour_patient(tour.sequence[2:end], patient_vec)
        end
        
    end

    sort!(c, by=x -> x.time, rev=true)
    outsiders = findall(x -> x == 0, counters[n_depot+1:end])
    for city in outsiders
        put_city_in_tour(c, city+n_depot, T, n_nodes, service)
    end

    #     chrm = Chromosome(Int[], 0.0, 0.0, c)
    #     for tour in c
    #         if tour.cost > chrm.fitness
    #             chrm.fitness = tour.cost
    #         end
    #         for city in tour.sequence
    #             push!(chrm.genes, city)
    #         end
    #     end
    #     return chrm

    child = Tour[]
    for tour in c
        if tour.time != 0.0 
            push!(child, tour)
        end
    end

    return child
end

function tour_crossover3(parent1::Chromosome, parent2::Chromosome, T::Matrix{Float64}, n_nodes::Int64, n_depot::Int64, distance::Matrix{Float64},
     service::Vector{Int64}, patient::Vector{Int64}, demand::Vector{Int64}, capacity_v::Int64, EC_G::Int64, VC_G::Int64)
    #3  At each step, select a tour from parent1, and select the tour with maximum mutual cities from parent2. 
    # Just keep the mutual cities in one of the tours
    # All the remaining cities will be placed in the current tours based on a greedy approach (minimum increase)
    P1 = deepcopy(parent1)
    P2 = deepcopy(parent2)
    c = Tour[]
    m = length(P1.tours)

    for i in 1:m
        tour1 = P1.tours[i].sequence[2:end]
        depot = tour1[1]
        max_intersection = -1
        tour2 = Int[]
        mutuals = Int[]
        r2 = 0
        for j in 1:length(P2.tours)
            mutual_cities = intersect(tour1, P2.tours[j].sequence[2:end])
            inter = length(mutual_cities)
            if inter > max_intersection
                max_intersection = inter
                tour2 = P2.tours[j].sequence[2:end]
                r2 = j
                mutuals = mutual_cities
            end
        end
        deleteat!(P2.tours, r2)
        mutual_indices = Int[]
        if rand() < 0.5
            for (j, node) in enumerate(tour1)
                if node in mutuals
                    push!(mutual_indices, j)
                end
            end
            cc = tour1[mutual_indices]
            staff = find_tour_staff(cc, demand)
            push!(c, Tour(cc, find_tour_length(cc, T, depot, service), staff, find_tour_patient(cc, patient), find_tour_operation_cost(cc, distance, staff, EC_G, VC_G)))
        else
            for (j, node) in enumerate(tour2)
                if node in mutuals
                    push!(mutual_indices, j)
                end
            end
            cc = tour2[mutual_indices]
            staff = find_tour_staff(cc, demand)
            push!(c, Tour(cc, find_tour_length(cc, T, depot, service), staff, find_tour_patient(cc, patient), find_tour_operation_cost(cc, distance, staff, EC_G, VC_G)))
        end
    end

    counters = zeros(n_nodes)
    outsiders = Int[]
    for tour in c
        delete_indices = Int[]
        for (j, node) in enumerate(tour.sequence[2:end])
            if counters[node] > 0
                push!(delete_indices, j)
            else
                counters[node] += 1
            end
        end
        if length(delete_indices) == length(tour.sequence[2:end]) # all of tour nodes are duplicated node 
            tour.time = 0.0
            tour.sequence = Int[tour.sequence[1]]
        elseif length(delete_indices) > 0 # subsection duplicated node 
            # tour.time = remove_cities_from_one_tour(tour, delete_indices, T)
            # deleteat!(tour.sequence, delete_indices.+1)
            # tour.time += sum([service[node] for node in tour.sequence[2:end]])
            deleteat!(tour.sequence, delete_indices.+1)
            tour.time = find_tour_length(tour.sequence[2:end], T, tour.sequence[1], service)
            tour.staff = find_tour_staff(tour.sequence, demand)
            tour.patient = find_tour_patient(tour.sequence[2:end], patient)
        end
    end

    sort!(c, by=x -> x.time, rev=true)
    outsiders = findall(x -> x == 0, counters[n_depot+1:end])
    for city in outsiders
        put_city_in_tour(c, city+n_depot, T, n_nodes)
    end

    #     chrm = Chromosome(Int[], 0.0, 0.0, c)
    #     for tour in c
    #         if tour.cost > chrm.fitness
    #             chrm.fitness = tour.cost
    #         end
    #         for city in tour.sequence
    #             push!(chrm.genes, city)
    #         end
    #     end
    #     return chrm

    child = Tour[]
    for tour in c
        if tour.time != 0.0 
            push!(child, tour)
        end
    end

    return child
end


function tour_crossover4(parent1::Chromosome, parent2::Chromosome, T::Matrix{Float64}, n_nodes::Int64, n_depot::Int64, distance::Matrix{Float64},
     service::Vector{Int64}, patient::Vector{Int64}, demand::Vector{Int64}, capacity_v::Int64, EC_G::Int64, VC_G::Int64)
    #10  At each step, select a tour from parent1, and select the tour with maximum mutual cities from parent2. 
    # Take the order of mutual cities from one tour and enforce it to the other one
    # All the remaining cities will be placed in the current tours based on a greedy approach (minimum increase)
    P1 = deepcopy(parent1)
    P2 = deepcopy(parent2)
    c = Tour[]
    m = length(P1.tours)

    for i in 1:m
        tour1 = P1.tours[i].sequence[2:end]
        depot = tour1[1]
        max_intersection = -1
        tour2 = Int[]
        mutuals = Int[]
        r2 = 0
        for j in 1:length(P2.tours)
            mutual_cities = intersect(tour1, P2.tours[j].sequence[2:end])
            inter = length(mutual_cities)
            if inter > max_intersection
                max_intersection = inter
                tour2 = P2.tours[j].sequence[2:end]
                r2 = j
                mutuals = mutual_cities
            end
        end
        deleteat!(P2.tours, r2)
        mutual_indices1 = Int[]
        mutual_indices2 = Int[]

        for (j, node) in enumerate(tour1)
            if node in mutuals
                push!(mutual_indices1, j)
            end
        end
        for (j, node) in enumerate(tour2)
            if node in mutuals
                push!(mutual_indices2, j)
            end
        end
        tour1[mutual_indices1] = tour2[mutual_indices2]
        staff = find_tour_staff(tour1, demand)
        push!(c, Tour(tour1, find_tour_length(tour1, T, depot, service), staff, find_tour_patient(tour1, patient), find_tour_operation_cost(tour1, distance, staff, EC_G, VC_G)))
    end
    counters = zeros(n_nodes)
    outsiders = Int[]
    for tour in c
        delete_indices = Int[]
        for (j, node) in enumerate(tour.sequence[2:end])
            if counters[node] > 0
                push!(delete_indices, j)
            else
                counters[node] += 1
            end
        end
        if length(delete_indices) == length(tour.sequence[2:end]) # all of tour nodes are duplicated node 
            tour.time = 0.0
            tour.sequence = Int[tour.sequence[1]]
        elseif length(delete_indices) > 0 # subsection duplicated node 
            # tour.time = remove_cities_from_one_tour(tour, delete_indices, T)
            # deleteat!(tour.sequence, delete_indices.+1)
            # tour.time += sum([service[node] for node in tour.sequence[2:end]])
            deleteat!(tour.sequence, delete_indices.+1)
            tour.time = find_tour_length(tour.sequence[2:end], T, tour.sequence[1], service)
            tour.staff = find_tour_staff(tour.sequence, demand)
            tour.patient = find_tour_patient(tour.sequence[2:end], patient)
        end
    end
    outsiders = findall(x -> x == 0, counters[n_depot+1:end])
    for city in outsiders
        put_city_in_tour(c, city+n_depot, T, n_nodes)
    end

    child = Tour[]
    sort!(c, by=x -> x.time, rev=true)
    for tour in c
        if tour.time != 0.0 
            push!(child, tour)
        end
    end

    return child
end


function tour_crossover5(parent1::Chromosome, parent2::Chromosome, T::Matrix{Float64}, n_nodes::Int64, n_depot::Int64)
    #10  At each step, select a tour from parent1, and select the tour with maximum mutual cities from parent2. 
    # keep the mutual cities in one of the tours, combine the remaining ones and randomly choose half of them and place them greedily
    # All the remaining cities will be placed in the current tours based on a greedy approach (minimum increase)
    P1 = deepcopy(parent1)
    P2 = deepcopy(parent2)
    c = Tour[]
    m = length(P1.tours)

    for i in 1:m
        tour1 = P1.tours[i].sequence
        depot = tour1[1]
        max_intersection = -1
        tour2 = Int[]
        mutuals = Int[]
        r2 = 0
        for j in 1:length(P2.tours)
            mutual_cities = intersect(tour1, P2.tours[j].sequence)
            inter = length(mutual_cities)
            if inter > max_intersection
                max_intersection = inter
                tour2 = P2.tours[j].sequence
                r2 = j
                mutuals = mutual_cities
            end
        end
        deleteat!(P2.tours, r2)
        mutual_indices1 = Int[]
        mutual_indices2 = Int[]
        rest = union(setdiff(tour1, mutuals), setdiff(tour2, mutuals))

        if length(rest) == 0
            if rand() < 0.5
                push!(c, Tour(tour1, find_tour_length(tour1, T, depot)))
            else
                push!(c, Tour(tour2, find_tour_length(tour2, T, depot)))
            end
        else
            rest_ = rest[1:rand(1:length(rest))]
            if rand() < 0.5
                for (j, node) in enumerate(tour1)
                    if node in mutuals
                        push!(mutual_indices1, j)
                    end
                end
                cc = copy(tour1[mutual_indices1])
            else
                for (j, node) in enumerate(tour2)
                    if node in mutuals
                        push!(mutual_indices2, j)
                    end
                end
                cc = copy(tour2[mutual_indices2])
            end
            ccc = Tour[]
            push!(ccc, Tour(cc, find_tour_length(cc, T, depot)))
            for city in rest_
                put_city_in_tour(ccc, city, T, n_nodes)
            end            
            push!(c, ccc[1])
        end
    end

    counters = zeros(n_nodes)
    outsiders = Int[]
    for tour in c
        delete_indices = Int[]
        for (j, node) in enumerate(tour.sequence[2:end])
            if counters[node] > 0
                push!(delete_indices, j)
            else
                counters[node] += 1
            end
        end
        if length(delete_indices) == length(tour.sequence[2:end]) # all of tour nodes are duplicated node 
            tour.cost = 0.0
            tour.sequence = Int[tour.sequence[1]]
        elseif length(delete_indices) > 0 # subsection duplicated node 
            tour.cost = remove_cities_from_one_tour(tour, delete_indices, T)
            deleteat!(tour.sequence, delete_indices.+1)
            tour.time += sum([service[node] for node in tour.sequence[2:end]])
        end
    end
    sort!(c, by=x -> x.cost, rev=true)
    outsiders = findall(x -> x == 0, counters[n_depot+1:end])
    for city in outsiders
        put_city_in_tour(c, city+n_depot, T, n_nodes)
    end

    child = Tour[]
    for tour in c
        if tour.cost != 0.0 
            push!(child, tour)
        end
    end

    return child
end

