function segments_intersect_orient(A::AbstractVector{Float64}, B::AbstractVector{Float64}, C::AbstractVector{Float64}, D::AbstractVector{Float64})
    # Unpack the coordinates of the first segment into separate variables
    x1, y1 = A
    x2, y2 = B

    # Unpack the coordinates of the second segment into separate variables
    x3, y3 = C
    x4, y4 = D

    # Calculate the orientations of the line segments
    o1 = orientation((x1, y1), (x2, y2), (x3, y3))
    o2 = orientation((x1, y1), (x2, y2), (x4, y4))
    o3 = orientation((x3, y3), (x4, y4), (x1, y1))
    o4 = orientation((x3, y3), (x4, y4), (x2, y2))

    # Check if the line segments intersect
    if o1 != o2 && o3 != o4
        return true
    end

    # Check if the line segments are collinear and overlapping
    #     if o1 == 0 && on_segment((x1, y1), (x3, y3), (x2, y2))
    #         return true
    #     end

    #     if o2 == 0 && on_segment((x1, y1), (x4, y4), (x2, y2))
    #         return true
    #     end

    #     if o3 == 0 && on_segment((x3, y3), (x1, y1), (x4, y4))
    #         return true
    #     end

    #     if o4 == 0 && on_segment((x3, y3), (x2, y2), (x4, y4))
    #         return true
    #     end

    # If none of the above conditions are met, the line segments do not intersect
    return false
end

function orientation(p, q, r)
    # Calculate the orientation of the triplet (p, q, r)
    val = (q[2] - p[2]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[2] - q[2])

    if val == 0
        # Points are collinear
        return 0
    elseif val > 0
        # Clockwise orientation
        return 1
    else
        # Counterclockwise orientation
        return 2
    end
end

function on_segment(p, q, r)
    # Check if point q lies on segment pr
    if q[0] <= max(p[0], r[0]) && q[0] >= min(p[0], r[0]) &&
       q[1] <= max(p[1], r[1]) && q[1] >= min(p[1], r[1])
        return true
    end

    return false
end


function find_intersections(t1::Vector{Int}, d1::Int64, t2::Vector{Int}, d2::Int64, Customers::Matrix{Float64}, Depots::Matrix{Float64})
    #     intersection_points = Vector{Vector{Int}}()
    
    nd = size(Depots)[1]
    for i = 0:length(t1)
        if i == 0 || i == length(t1)
            loopset = [k for k in 1:length(t2)-1]
        else
            loopset = [k for k = 0:length(t2)]
        end

        if i == 0
            node1 = @view Depots[d1, :]  # essentially same as `node1 = depot`, but for type stability
            node2 = @view Customers[t1[i+1]-nd, :]
        elseif i == length(t1)
            node1 = @view Customers[t1[i]-nd, :]
            node2 = @view Depots[d1, :]
        else
            node1 = @view Customers[t1[i]-nd, :]
            node2 = @view Customers[t1[i+1]-nd, :]
        end    

        for j in loopset
            intersected = false

            if j == 0
                node3 = @view Depots[d2, :]
                node4 = @view Customers[t2[j+1]-nd, :]
            elseif j == length(t2)
                node3 = @view Customers[t2[j]-nd, :]
                node4 = @view Depots[d2, :]
            else
                node3 = @view Customers[t2[j]-nd, :]
                node4 = @view Customers[t2[j+1]-nd, :]
            end

            if segments_intersect_orient(node1, node2, node3, node4)
                # println(i," ",j)
                # println(node1," ", node2," ",node3," ",node4)
                return i, j, true
            end
        end
    end
    return 0, 0, false
end

function solve_one_intersections(t1::Vector{Int}, d1::Int64, t2::Vector{Int}, d2::Int64, T::Matrix{Float64}, k1::Int, k2::Int, service::Vector{Int64})
    n1 = length(t1)
    n2 = length(t2)
    if k1 == 0
        a = Int[]
    else
        a = t1[1:k1]
    end
    if k2 == n2
        b = Int[]
    else
        b = t2[k2+1:n2]
    end
    if k2 == 0
        c = Int[]
    else
        c = t2[1:k2]
    end
    if k1 == n1
        d = Int[]
    else
        d = t1[k1+1:n1]
    end
    # Will be adding operation, which is change depots each tour 
    tour11 = vcat(a, b)
    tour12 = vcat(c, d)
    tour21 = vcat(a, reverse(c))
    tour22 = vcat(reverse(b), d)
    time11 = find_tour_length(tour11, T, d1, service)
    time12 = find_tour_length(tour12, T, d2, service)
    time21 = find_tour_length(tour21, T, d1, service)
    time22 = find_tour_length(tour22, T, d2, service)
    if max(time11, time12) < max(time21, time22)
        return tour11, tour12, time11, time12
    else
        return tour21, tour22, time21, time22
    end
end

function solve_all_intersections!(chrm::Chromosome, Customers::Matrix{Float64}, Depots::Matrix{Float64}, T::Matrix{Float64}, 
    service::Vector{Int64}, patient::Vector{Int64}, demand::Vector{Int64}, capacity_v::Int64)
    # println("Start solve_all_intersections!")
    m = length(chrm.tours)
    n_nodes = length(chrm.genes)
    intersected_ = true
    cnt = 0
    while intersected_
        intersected_ = false
        for i in 1:m-1
            for j = i+1:m
                tour1 = chrm.tours[i].sequence[2:end]
                depot1 = chrm.tours[i].sequence[1]
                tour2 = chrm.tours[j].sequence[2:end]
                depot2 = chrm.tours[j].sequence[1]
                if length(tour1) > 0 && length(tour2) > 0
                    intersected = true
                    while intersected
                        k1, k2, intersected = find_intersections(tour1, depot1, tour2, depot2, Customers, Depots)
                        if cnt >= 100 intersected = false end
                        if intersected
                            cnt += 1
                            t1, t2, c1, c2 = solve_one_intersections(tour1, depot1, tour2, depot2, T, k1, k2, service)

                            if capacity_v - find_tour_staff(t1, demand) < sum(find_tour_patient(t1, patient)) break end
                            if capacity_v - find_tour_staff(t2, demand) < sum(find_tour_patient(t2, patient)) break end

                            tour1 = copy(t1)
                            tour2 = copy(t2)
                            chrm.tours[i].time = c1
                            chrm.tours[i].sequence = vcat(depot1, tour1)
                            chrm.tours[i].staff = find_tour_staff(chrm.tours[i].sequence, demand)
                            chrm.tours[i].patient = find_tour_patient(chrm.tours[i].sequence[2:end], patient)
                            # chrm.tours[i].operation_cost = find_tour_operation_cost(chrm.tours[i].sequence, distance, staff, EC_G, VC_G)

                            chrm.tours[j].time = c2
                            chrm.tours[j].sequence = vcat(depot2, tour2)
                            chrm.tours[j].staff = find_tour_staff(chrm.tours[j].sequence, demand)
                            chrm.tours[j].patient = find_tour_patient(chrm.tours[j].sequence[2:end], patient)
                            # chrm.tours[j].operation_cost = find_tour_operation_cost(chrm.tours[j].sequence, distance, staff, EC_G, VC_G)
                            intersected_ = true
                        end
                    end
                end
            end
        end
        for tour in chrm.tours
            two_opt_on_route(tour, T, n_nodes)
        end
    end

    
    if rand() < 0.1
        improve_after_removing_intersections(chrm.tours, T, n_nodes, m, Customers, Depots, service, patient, demand, capacity_v)
    end
    chrm.genes = Int[]
    chrm.fitness = 0.0
    for tour in chrm.tours
        chrm.genes = vcat(chrm.genes, tour.sequence, tour.sequence[1], -1)
        if chrm.fitness < tour.time
            chrm.fitness = tour.time
        end
        # tour.operation_cost = 0.0
    end

    # chrm.total_cost = calculate_costs(chrm.tours, distance, FC, EC_G, VC_G) + r_cost
end

function two_opt_on_route(tour::Tour, T::Matrix{Float64}, n_nodes::Int)   #2-opt 
    tour1 = tour.sequence[2:end]
    depot1 = tour.sequence[1]
    time1 = tour.time
    nt = length(tour1)
    improved = true
    while improved
        improved = false
        for i1 in 1:nt-1
            if improved
                break
            end
            for i2 = i1+1:nt
                new_time = calculate_new_cost_2_opt(tour1, depot1, time1, i1, i2, T, n_nodes)
                if round(new_time, digits=2) < round(time1, digits=2)
                    time1 = new_time
                    improved = true
                    tour1[i1:i2] = reverse(tour1[i1:i2])
                    tour.time = new_time
                    tour.sequence = vcat(depot1, tour1)
                    break
                end
            end
        end
    end
end

function improve_after_removing_intersections(tours::Vector{Tour}, T::Matrix{Float64}, n_nodes::Int, m::Int, Customers::Matrix{Float64}, Depots::Matrix{Float64},
    service::Vector{Int64}, patient::Vector{Int64}, demand::Vector{Int64}, capacity_v::Int64)

    nd = size(Depots)[1]
    
    sort!(tours, by=x -> x.time, rev=true)
    means = [mean(Customers[t1[2:end].-nd, :], dims=1)[1, :] for t1 in [tours[i].sequence for i in 1:m]]

    distances = zeros(m, m)
    for i in 1:m-1
        for j = i+1:m
            distances[i, j] = euclidean(means[i], means[j])
            distances[j, i] = distances[i, j]
        end
    end

    all_tours = [i for i in 1:m]

    improved = true

    while improved
        current_tour = 1
        improved = false
        next_tours = sortperm(distances[current_tour, :])[2:min(m, 3)]
        #         println(tours[current_tour].cost)
        for next_tour in next_tours
            #             println(tours[next_tour].cost)
            tour1 = tours[current_tour].sequence[2:end]
            depot1 = tours[current_tour].sequence[1]
            time1 = tours[current_tour].time

            tour2 = tours[next_tour].sequence[2:end]
            depot2 = tours[next_tour].sequence[1]
            time2 = tours[next_tour].time
            for i in 1:length(tour1)
                # Capa check
                if capacity_v - max(tours[next_tour].staff, demand[tour1[i]]) < sum(tours[next_tour].patient) + patient[tour1[i]] continue end
                new_time1 = calculate_new_cost_remove_one(tour1, depot1, time1, i, T) - service[tour1[i]]
                new_time2 = Inf
                best_position = 0
                for j in 1:length(tour2)+1
                    temp = calculate_new_cost_add_one(tour2, depot2, time2, tour1[i], j, T) + service[tour1[i]]
                    if temp < new_time2
                        new_time2 = temp
                        best_position = j
                    end
                end
                #                 println(i,"  ", best_position)

                if new_time2 < time1
                    t2 = copy(tour2)
                    insert!(t2, best_position, tour1[i])
                    if i == 1
                        node1 = Depots[depot1, :]
                        node2 = Customers[tour1[i]-nd, :] # node2 = Customers[tour1[i+1]-nd, :] ==> node2 = Customers[tour1[i]-nd, :]. 24-07-11. JH 
                    elseif i == length(tour1)
                        node1 = Customers[tour1[i-1]-nd, :]
                        node2 = Depots[depot1, :]
                    else
                        node1 = Customers[tour1[i-1]-nd, :]
                        node2 = Customers[tour1[i+1]-nd, :]
                    end
                    node3 = Customers[tour1[i]-nd, :]
                    intersected = false
                    for j in 1:length(t2)
                        if j == 1
                            node4 = Depots[depot2, :]
                            node5 = Customers[t2[j]-nd, :]
                        elseif j == length(t2) + 1
                            node4 = Customers[t2[length(t2)]-nd, :]
                            node5 = Depots[depot2, :]
                        else
                            node4 = Customers[t2[j-1]-nd, :]
                            node5 = Customers[t2[j]-nd, :]
                        end

                        intersect1 = segments_intersect_orient(node1, node2, node3, node4)
                        intersect2 = segments_intersect_orient(node1, node2, node3, node5)
                        if intersect1 || intersect2
                            intersected = true
                            break
                        end
                    end
                    if !intersected
                        insert!(tour2, best_position, tour1[i])
                        deleteat!(tour1, i)
                        tours[current_tour].time = new_time1
                        tours[current_tour].sequence = vcat(depot1, tour1)
                        tours[current_tour].staff = find_tour_staff(tours[current_tour].sequence, demand)
                        tours[current_tour].patient = find_tour_patient(tours[current_tour].sequence[2:end], patient)
                        # tours[current_tour].operation_cost = find_tour_operation_cost(tours[current_tour].sequence, distance, staff, EC_G, VC_G)

                        tours[next_tour].time = new_time2
                        tours[next_tour].sequence = vcat(depot2, tour2)
                        tours[next_tour].staff = find_tour_staff(tours[next_tour].sequence, demand)
                        tours[next_tour].patient = find_tour_patient(tours[next_tour].sequence[2:end], patient)
                        # tours[next_tour].operation_cost = find_tour_operation_cost(tours[next_tour].sequence, distance, staff, EC_G, VC_G)

                        sort!(tours, by=x -> x.time, rev=true)
                        improved = true
                        break
                    end
                end
            end
        end
    end
end