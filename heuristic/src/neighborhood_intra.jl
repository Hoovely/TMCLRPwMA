# n_nodes+2 ==> depot, n_nodes+1 ?? (depot)
function Ni1!(chrm::Chromosome, TT::Matrix{Float64}, close_nodes::Matrix{Bool}, n_nodes::Int)   #Reinsert
    r1 = 1
    if rand() < 0.5
        r1 = argmax([chrm.tours[i].time for i in 1:length(chrm.tours)])
    else
        r1 = rand(1:length(chrm.tours))
    end

    tour1 = chrm.tours[r1].sequence[2:end]
    depot1 = chrm.tours[r1].sequence[1]
    time1 = chrm.tours[r1].time


    nt = length(tour1)
    if nt <= 1
        return
    end

    k1 = rand(1:length(tour1))
    city1 = tour1[k1]
    candidates = Int[]
    if nt == 2
        candidates = [1, 2]
    else
        if close_nodes[depot1, city1] || close_nodes[city1, tour1[1]]
            push!(candidates, 1)
        end
        for i = 2:nt-1
            if i > k1
                if close_nodes[city1, tour1[i]] || close_nodes[city1, tour1[i+1]]
                    push!(candidates, i)
                end
            elseif i < k1
                if close_nodes[city1, tour1[i]] || close_nodes[city1, tour1[i-1]]
                    push!(candidates, i)
                end
            end
        end
        if close_nodes[depot1, city1] || close_nodes[city1, tour1[nt]]
            push!(candidates, nt)
        end
    end
    candidates = collect(setdiff(Set(candidates), Set([k1])))
    if length(candidates) == 0
        return
    end

    k2 = candidates[rand(1:length(candidates))]
    #     k2 = rand(1:length(tour1))
    new_time1 = calculate_new_cost_exchange_one(tour1, depot1, time1, city1, k1, k2, TT, n_nodes)

    if new_time1 >= time1
        return
    end

    deleteat!(tour1, k1)
    insert!(tour1, k2, city1)

    chrm.tours[r1].time = new_time1
    chrm.tours[r1].sequence = vcat(depot1, tour1)
    chrm.fitness = maximum([chrm.tours[i].time for i in 1:length(chrm.tours)])
    index = 0
    for tour in chrm.tours
        chrm.genes[index+1:index+length(tour.sequence)+2] = vcat(tour.sequence, tour.sequence[1], -1)
        index += length(tour.sequence)+2
    end
end


function Ni2!(chrm::Chromosome, TT::Matrix{Float64}, close_nodes::Matrix{Bool}, service::Vector{Int})   #Exchange (permutation between two customers)
    r1 = 1
    if rand() < 0.5
        r1 = argmax([chrm.tours[i].time for i in 1:length(chrm.tours)])
    else
        r1 = rand(1:length(chrm.tours))
    end
    tour1 = chrm.tours[r1].sequence[2:end]
    depot1 = chrm.tours[r1].sequence[1]
    if length(tour1) <= 1
        return
    end
    time1 = chrm.tours[r1].time
    nt = length(tour1)
    k1 = rand(1:nt)
    city1 = tour1[k1]

    candidates = Int[]

    for i in 1:nt
        if i != k1
            if i == 1
                if close_nodes[depot1, city1] || close_nodes[city1, tour1[2]]
                    push!(candidates, 1)
                end
            elseif i == nt
                if close_nodes[depot1, city1] || close_nodes[city1, tour1[nt]]
                    push!(candidates, nt)
                end
            else
                if close_nodes[city1, tour1[i-1]] || close_nodes[city1, tour1[i+1]]
                    push!(candidates, i)
                end
            end
        end
    end
    if length(candidates) == 0
        return
    end
    k2 = candidates[rand(1:length(candidates))]
    #     k2 = rand(1:nt)
    city2 = tour1[k2]

    new_time1 = calculate_new_cost_exchange_two(tour1, depot1, time1, city1, k1, city2, k2, TT, service)

    if new_time1 >= time1
        return
    end

    tour1[k1] = city2
    tour1[k2] = city1
    chrm.tours[r1].time = new_time1
    chrm.tours[r1].sequence = vcat(depot1, tour1)
    chrm.fitness = maximum([chrm.tours[i].time for i in 1:length(chrm.tours)])
    index = 0
    for tour in chrm.tours
        chrm.genes[index+1:index+length(tour.sequence)+2] = vcat(tour.sequence, tour.sequence[1], -1)
        index += length(tour.sequence)+2
    end
end

function Ni3!(chrm::Chromosome, T::Matrix{Float64}, close_nodes::Matrix{Bool}, n_nodes::Int)   #Or-opt2 
    r1 = 1
    if rand() < 0.5
        r1 = argmax([chrm.tours[i].time for i in 1:length(chrm.tours)])
    else
        r1 = rand(1:length(chrm.tours))
    end
    tour1 = chrm.tours[r1].sequence[2:end]
    depot1 = chrm.tours[r1].sequence[1]
    if length(tour1) <= 2
        return
    end
    time1 = chrm.tours[r1].time
    nt = length(tour1)
    k1 = rand(1:nt-1)
    city1 = tour1[k1]
    city2 = tour1[k1+1]

    candidates = Int[]
    if nt == 3
        candidates = [1, 2]
    else
        if close_nodes[depot1, city1] || close_nodes[city2, tour1[1]]
            push!(candidates, 1)
        end
        for i = 2:nt-2
            if i > k1
                if close_nodes[city1, tour1[i+1]] || close_nodes[city2, tour1[i+2]]
                    push!(candidates, i)
                end
            elseif i < k1
                if close_nodes[city1, tour1[i-1]] || close_nodes[city2, tour1[i]]
                    push!(candidates, i)
                end
            end
        end
        if close_nodes[depot1, city2] || close_nodes[city1, tour1[nt]]
            push!(candidates, nt - 1)
        end
    end
    candidates = collect(setdiff(Set(candidates), Set([k1])))
    if length(candidates) == 0
        return
    end

    #     k2 = rand(1:length(tour1)-1)   #Way to improve 
    k2 = candidates[rand(1:length(candidates))]

    z1 = calculate_new_cost_or_opt2(tour1, depot1, time1, city1, k1, city2, k2, T, n_nodes)

    if z1 >= time1
        return
    end
    deleteat!(tour1, [k1, k1 + 1])
    insert!(tour1, k2, city1)
    insert!(tour1, k2 + 1, city2)
    chrm.tours[r1].time = z1
    chrm.tours[r1].sequence = vcat(depot1, tour1)
    chrm.fitness = maximum([chrm.tours[i].time for i in 1:length(chrm.tours)])
    index = 0
    for tour in chrm.tours
        chrm.genes[index+1:index+length(tour.sequence)+2] = vcat(tour.sequence, tour.sequence[1], -1)
        index += length(tour.sequence)+2
    end

    return
end


function Ni4!(chrm::Chromosome, T::Matrix{Float64}, close_nodes::Matrix{Bool}, n_nodes::Int)   #Or-opt3 
    r1 = 1
    if rand() < 0.5
        r1 = argmax([chrm.tours[i].time for i in 1:length(chrm.tours)])
    else
        r1 = rand(1:length(chrm.tours))
    end
    tour1 = chrm.tours[r1].sequence[2:end]
    depot1 = chrm.tours[r1].sequence[1]
    if length(tour1) <= 3
        return
    end
    time1 = chrm.tours[r1].time
    nt = length(tour1)
    k1 = rand(1:nt-2)
    city1 = tour1[k1]
    city2 = tour1[k1+1]
    city3 = tour1[k1+2]
    candidates = Int[]
    if nt == 4
        candidates = [1, 2]
    else
        if close_nodes[depot1, city1] || close_nodes[city3, tour1[1]]
            push!(candidates, 1)
        end
        for i = 2:nt-3
            if i > k1
                if close_nodes[city1, tour1[i+2]] || close_nodes[city3, tour1[i+3]]
                    push!(candidates, i)
                end
            elseif i < k1
                if close_nodes[city1, tour1[i-1]] || close_nodes[city3, tour1[i]]
                    push!(candidates, i)
                end
            end
        end
        if close_nodes[depot1, city3] || close_nodes[city1, tour1[nt]]
            push!(candidates, nt - 2)
        end
    end
    candidates = collect(setdiff(Set(candidates), Set([k1])))
    if length(candidates) == 0
        return
    end

    k2 = candidates[rand(1:length(candidates))]
    new_time1 = calculate_new_cost_or_opt3(tour1, depot1, time1, city1, city2, city3, k1, k2, T, n_nodes)

    if new_time1 >= time1
        return
    end
    deleteat!(tour1, [k1, k1 + 1, k1 + 2])
    insert!(tour1, k2, city1)
    insert!(tour1, k2 + 1, city2)
    insert!(tour1, k2 + 2, city3)
    chrm.tours[r1].time = new_time1
    chrm.tours[r1].sequence = vcat(depot1, tour1)
    chrm.fitness = maximum([chrm.tours[i].time for i in 1:length(chrm.tours)])
    index = 0
    for tour in chrm.tours
        chrm.genes[index+1:index+length(tour.sequence)+2] = vcat(tour.sequence, tour.sequence[1], -1)
        index += length(tour.sequence)+2
    end

    return
end

function Ni5!(chrm::Chromosome, T::Matrix{Float64}, close_nodes::Matrix{Bool}, n_nodes::Int)   #2-opt 
    r1 = 1
    if rand() < 0.5
        r1 = argmax([chrm.tours[i].time for i in 1:length(chrm.tours)])
    else
        r1 = rand(1:length(chrm.tours))
    end
    tour1 = chrm.tours[r1].sequence[2:end]
    depot1 = chrm.tours[r1].sequence[1]
    if length(tour1) <= 2
        return
    end
    time1 = chrm.tours[r1].time
    nt = length(tour1)
    # indices = sample(1:length(tour1), 2, replace = false)
    i1 = rand(1:length(tour1))

    candidates = Int[]

    if i1 == 1
        for i = 2:nt-1
            if close_nodes[depot1, tour1[i]] || close_nodes[tour1[i1+1], tour1[1]]
                push!(candidates, i)
            end
        end
    elseif i1 == nt
        for i = 2:nt-1
            if close_nodes[depot1, tour1[i]] || close_nodes[tour1[i1-1], tour1[nt]]
                push!(candidates, i)
            end
        end
    else
        if close_nodes[depot1, tour1[i1]] || close_nodes[tour1[i1+1], tour1[1] ]
            push!(candidates, 1)
        end
        if close_nodes[depot1, tour1[i1]] || close_nodes[tour1[i1-1], tour1[nt]]
            push!(candidates, nt)
        end
        for i = 2:nt-1
            if i > i1
                if close_nodes[tour1[i+1], tour1[i1]] || close_nodes[tour1[i1-1], tour1[i]]
                    push!(candidates, i)
                end
            elseif i < i1
                if close_nodes[tour1[i-1], tour1[i1]] || close_nodes[tour1[i1+1], tour1[i]]
                    push!(candidates, i)
                end
            end
        end
    end

    if length(candidates) == 0
        return
    end
    candidates = collect(Set(candidates))
    i2 = candidates[rand(1:length(candidates))]
    k1, k2 = min(i1, i2), max(i1, i2)
    new_time = calculate_new_cost_2_opt(tour1, depot1, time1, k1, k2, T, n_nodes)

    if new_time >= time1
        return
    end
    tour1[k1:k2] = reverse(tour1[k1:k2])
    chrm.tours[r1].time = new_time
    chrm.tours[r1].sequence = vcat(depot1, tour1)
    chrm.fitness = maximum([chrm.tours[i].time for i in 1:length(chrm.tours)])
    index = 0
    for tour in chrm.tours
        chrm.genes[index+1:index+length(tour.sequence)+2] = vcat(tour.sequence, tour.sequence[1], -1)
        index += length(tour.sequence)+2
    end

    return
end

function Ni6!(chrm::Chromosome, T::Matrix{Float64}, close_nodes::Matrix{Bool}, n_nodes::Int)   #3-opt 
    r1 = 1
    if rand() < 0.5
        r1 = argmax([chrm.tours[i].time for i in 1:length(chrm.tours)])
    else
        r1 = rand(1:length(chrm.tours))
    end
    tour1 = chrm.tours[r1].sequence[2:end]
    depot1 = chrm.tours[r1].sequence[1]
    if length(tour1) <= 2
        return
    end
    time1 = chrm.tours[r1].time
    nt = length(tour1)
    k1, k2, k3 = sort!(sample(1:length(tour1), 3, replace=false))

    new_time = calculate_new_cost_3_opt(tour1, depot1, time1, k1, k2, k3, T, n_nodes)
    if new_time >= time1
        return
    end
    if k2 - k1 >= 3
        tour1[k1+1:k2-1] = reverse(tour1[k1+1:k2-1])
    end
    if k3 - k2 >= 3
        tour1[k2+1:k3-1] = reverse(tour1[k2+1:k3-1])
    end
    chrm.tours[r1].time = new_time
    chrm.tours[r1].sequence = vcat(depot1, tour1)
    chrm.fitness = maximum([chrm.tours[i].time for i in 1:length(chrm.tours)])
    index = 0
    for tour in chrm.tours
        l = length(tour.sequence)
        chrm.genes[index+1:index+l] = tour.sequence
        index += l
    end
end

function Ni7!(chrm::Chromosome, T::Matrix{Float64}, close_nodes::Matrix{Bool}, n_nodes::Int)   #3-permute 
    r1 = 1
    if rand() < 0.5
        r1 = argmax([chrm.tours[i].time for i in 1:length(chrm.tours)])
    else
        r1 = rand(1:length(chrm.tours))
    end
    tour1 = chrm.tours[r1].sequence
    if length(tour1) <= 2
        return
    end
    time1 = chrm.tours[r1].time
    nt = length(tour1)
    k1 = rand(1:length(tour1)-2)
    temp1 = copy(tour1[k1:k1+2])
    temp2 = shuffle(temp1)

    new_time = calculate_new_cost_3_permute(tour1, time1, temp1, temp2, k1, T, n_nodes)

    if new_time >= time1
        return
    end
    tour1[k1:k1+2] = temp2
    chrm.tours[r1].time = new_time
    chrm.fitness = maximum([chrm.tours[i].time for i in 1:length(chrm.tours)])
    index = 0
    for tour in chrm.tours
        l = length(tour.sequence)
        chrm.genes[index+1:index+l] = tour.sequence
        index += l
    end

end