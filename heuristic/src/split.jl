mutable struct Label
    Ri::Vector{Int}
    Vir::Vector{Float64}
    Pir::Vector{Int}
    Cir::Vector{Float64}
end


function Get_feasible_K(S::Vector{Int64}, Q::Int64, patient::Vector{Int64}, demand::Vector{Int64})

    # K = ceil(Int64, sum([patient[node] for node in S])/(Q-maximum([demand[node] for node in S]))) + 3
    K = ceil(Int64, length(S)//2)

    return K
end

"
    SPLIT func
    Description: 
        - Divide MDTSP by MDVRP
    Params:
        - Travel time matrix
        - number of vehicle
        - MDTSP tours vector 
        - number of node
        - service time vector
        - number of patient vector
        - number of demand vector
    Return:
        - Vector of red patient allocated hospital, Rescue completion time
"
function SPLIT(TT::Matrix{Float64}, tours::Vector{Vector{Int64}}, service::Vector{Int64}, patient::Vector{Int64}, demand::Vector{Int64}, capacity_v::Int64)

    obj_vec = Vector{Float64}()
    trips_vec = Vector{Tour}()

    #In m-TSP, demands is a vector of ones and W is infinity
    for tour in tours
        # println("tour is $tour")
        depot = tour[1]
        S = tour[2:end]
        K = Get_feasible_K(S, capacity_v, patient, demand)
        # println(K)
        n = length(S)
        labels = Vector{Label}(undef, n)

        @inbounds for i in 1:n
            R = Int[]
            if i == n
                R = [j for j in 1:K]
            else
                R = [j for j in 1:min(i, K - 1)]
            end
            V = fill(Inf, length(R))
            P = fill(n+1, length(R))
            C = fill(Inf, length(R))
            labels[i] = Label(R, V, P, C)
        end

        @inbounds for i in 1:n
            R = Int[0]
            if i > 1
                R = labels[i-1].Ri
            end
            for r in R
                Current_V = 0.0
                if i > 1
                    Current_V = labels[i].Vir[r]
                end
                if Current_V < Inf    #This is V^i_r
                    load = 0
                    capa = capacity_v
                    t = zero(eltype(TT))
                    j = i
                    temp = 0
                    while ( (j <= n) && (load + patient[S[j]] <= capa + temp - max(temp, demand[S[j]])) )

                        load += patient[S[j]]
                        if i == j
                            t = TT[depot, S[j]] + service[S[j]]
                            capa -= demand[S[j]]
                            temp = demand[S[j]]
                        else
                            t = t + TT[S[j-1], S[j]] + service[S[j]]
                            capa += temp 
                            temp = max(temp, demand[S[j]])
                            capa -= temp 
                        end
                   
                        if r + 1 in labels[j].Ri
                            old_t = 0.0
                            if i > 1
                                old_t = labels[i-1].Vir[r]
                            end
                            new_t = max(old_t, t)
                            if new_t < labels[j].Vir[r+1]
                                labels[j].Vir[r+1] = new_t
                                labels[j].Pir[r+1] = i - 1
                                labels[j].Cir[r+1] = t
                            end
                        end

                        j += 1
                    end
                end
            end
        end
        #     return labels
        #     rs = argmin(labels[n].Vir)
        rs = K
        trips = Vector{Tour}(undef, rs)
        @inbounds for i in 1:rs
            trips[i] = Tour(Int[depot], 0.0, 0, Vector{Int}(), 0.0)
        end
        tt = rs
        j = n
        while tt > 0
            i = labels[j].Pir[tt]
            trips[tt].time = labels[j].Cir[tt]
            for k = i+1:j
                node = S[k]
                push!(trips[tt].sequence, node)
                trips[tt].staff = max(trips[tt].staff, demand[node])
                push!(trips[tt].patient, patient[node])
            end
            tt -= 1
            j = i
        end
        obj = minimum(labels[n].Vir)
        push!(obj_vec, obj)
        trips_vec = vcat(trips_vec, trips)
    end

    return maximum(obj_vec), trips_vec
end


mutable struct Label_
    Vir::Vector{Float64}
    Pir::Vector{Int}
    Cir::Vector{Float64}
end


function SPLIT_test(TT::Matrix{Float64}, K::Int, S::Vector{Int}) #In m-TSP, demands is a vector of ones and W is infinity
    n = size(TT)[1] - 2
    labels = Label_[]
    for i in 1:n
        V = Float64[]
        P = Int[]
        C = Float64[]
        for r in 1:K
            push!(V, Inf)
            push!(P, n + 1)
            push!(C, Inf)
        end
        push!(labels, Label_(V, P, C))
    end
    i = 1
    r = 0
    t = 0
    j = i
    while (j <= n)
        if i == j
            t = TT[1, S[j]+1] + TT[S[j]+1, 1]
        else
            t = t - TT[S[j-1]+1, 1] + TT[S[j-1]+1, S[j]+1] + TT[S[j]+1, 1]
        end
        old_t = 0.0
        new_t = max(old_t, t)
        if new_t < labels[j].Vir[r+1]
            labels[j].Vir[r+1] = new_t
            labels[j].Pir[r+1] = i - 1
            labels[j].Cir[r+1] = t
        end
        j += 1
    end

    for i = 2:n
        for r in 1:K-1
            Current_V = labels[i].Vir[r]
            if Current_V < Inf    #This is V^i_r
                t = 0
                j = i
                while (j <= n)
                    if i == j
                        t = TT[1, S[j]+1] + TT[S[j]+1, 1]
                    else
                        t = t - TT[S[j-1]+1, 1] + TT[S[j-1]+1, S[j]+1] + TT[S[j]+1, 1]
                    end
                    old_t = 0.0
                    if i > 1
                        old_t = labels[i-1].Vir[r]
                    end
                    new_t = max(old_t, t)
                    if new_t < labels[j].Vir[r+1]
                        labels[j].Vir[r+1] = new_t
                        labels[j].Pir[r+1] = i - 1
                        labels[j].Cir[r+1] = t
                    end
                    j += 1
                end
            end
        end
    end
    #     return labels

    trips = Vector{Tour}()

    rs = argmin(labels[n].Vir)
    #     rs = K
    for i in 1:rs
        push!(trips, Tour(Int[], 0.0))
    end
    t = rs
    j = n
    while t > 0
        i = labels[j].Pir[t]
        trips[t].cost = labels[j].Cir[t]
        for k = i+1:j
            push!(trips[t].sequence, S[k])
        end
        t -= 1
        j = i
    end
    obj = minimum(labels[n].Vir)
    return obj, trips
end