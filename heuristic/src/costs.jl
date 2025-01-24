function calculate_new_cost_add_one(tour::Vector{Int}, depot::Int64, cost::Float64, city::Int, position::Int, T::Matrix{Float64})
    nt = length(tour)
    if nt == 0
        return T[depot, city]
    else
        if position == 1
            cost += T[depot, city] + T[city, tour[1]] - T[depot, tour[1]]
        elseif position == nt + 1
            cost += T[tour[nt], city]
        else
            cost += T[tour[position-1], city] + T[city, tour[position]] - T[tour[position-1], tour[position]]
        end
    end
    return cost
end

function calculate_new_cost_add_two(tour::Vector{Int}, depot::Int64, cost::Float64, city1::Int, city2::Int, position::Int, T::Matrix{Float64}, n_nodes::Int)
    nt = length(tour)

    cost1 = cost
    cost2 = cost
    if nt == 0
        cost1 = T[depot, city1] + T[city1, city2] + T[city2, depot]
        cost2 = T[depot, city2] + T[city2, city1] + T[city1, depot]
    else
        if position == 1
            cost1 += T[depot, city1] + T[city1, city2] + T[city2, tour[1]] - T[depot, tour[1]]
            cost2 += T[depot, city2] + T[city2, city1] + T[city1, tour[1]] - T[depot, tour[1]]
        elseif position == nt + 1
            cost1 += T[city2, depot] + T[city1, city2] + T[tour[nt], city1] - T[tour[nt], depot]
            cost2 += T[city1, depot] + T[city2, city1] + T[tour[nt], city2] - T[tour[nt], depot]
        else
            cost1 += T[tour[position-1], city1] + T[city1, city2] + T[city2, tour[position]] - T[tour[position-1], tour[position]]
            cost2 += T[tour[position-1], city2] + T[city2, city1] + T[city1, tour[position]] - T[tour[position-1], tour[position]]
        end
    end
    if cost1 < cost2
        return cost1, true
    else
        return cost2, false
    end
end

function calculate_new_cost_add_three(tour::Vector{Int}, cost::Float64, city1::Int, city2::Int, city3::Int, position::Int, T::Matrix{Float64}, n_nodes::Int)
    nt = length(tour)
    cost1 = cost
    cost2 = cost
    if nt == 0
        cost1 = T[1, city1+1] + T[city1+1, city2+1] + T[city2+1, city3+1] + T[city3+1, n_nodes+2]
        cost2 = T[1, city3+1] + T[city3+1, city2+1] + T[city2+1, city1+1] + T[city3+1, n_nodes+2]
    else
        if position == 1
            cost1 += T[1, city1+1] + T[city1+1, city2+1] + T[city2+1, city3+1] + T[city3+1, tour[1]+1] - T[1, tour[1]+1]
            cost2 += T[1, city3+1] + T[city3+1, city2+1] + T[city2+1, city1+1] + T[city1+1, tour[1]+1] - T[1, tour[1]+1]
        elseif position == nt + 1
            cost1 += T[city3+1, n_nodes+2] + T[city1+1, city2+1] + T[city2+1, city3+1] + T[tour[nt]+1, city1+1] - T[tour[nt]+1, n_nodes+2]
            cost2 += T[city1+1, n_nodes+2] + T[city3+1, city2+1] + T[city2+1, city1+1] + T[tour[nt]+1, city3+1] - T[tour[nt]+1, n_nodes+2]
        else
            cost1 += T[tour[position-1]+1, city1+1] + T[city1+1, city2+1] + T[city2+1, city3+1] + T[city3+1, tour[position]+1] - T[tour[position-1]+1, tour[position]+1]
            cost2 += T[tour[position-1]+1, city3+1] + T[city3+1, city2+1] + T[city2+1, city1+1] + T[city1+1, tour[position]+1] - T[tour[position-1]+1, tour[position]+1]
        end
    end
    if cost1 < cost2
        return cost1, true
    else
        return cost2, false
    end
end

function calculate_new_cost_remove_one(tour::Vector{Int}, depot::Int64, cost::Float64, position::Int, T::Matrix{Float64})
    nt = length(tour)

    if nt == 1
        return 0.0
    else
        if position == 1
            cost += T[depot, tour[2]] - T[tour[1], tour[2]] - T[depot, tour[1]]
        elseif position == nt
            cost -= T[tour[nt-1], tour[nt]]
        else
            cost += T[tour[position-1], tour[position+1]] - T[tour[position-1], tour[position]] - T[tour[position], tour[position+1]]
        end
    end
    return cost
end

function calculate_new_cost_remove_two(tour::Vector{Int}, cost::Float64, position::Int, T::Matrix{Float64}, n_nodes::Int)
    nt = length(tour)
    if nt == 2
        return 0.0
    else
        if position == 1
            cost += T[1, tour[3]+1] - T[tour[1]+1, tour[2]+1] - T[tour[2]+1, tour[3]+1] - T[1, tour[1]+1]
        elseif position == nt - 1
            cost += T[tour[nt-2]+1, n_nodes+2] - T[tour[nt-2]+1, tour[nt-1]+1] - T[tour[nt-1]+1, tour[nt]+1] - T[tour[nt]+1, n_nodes+2]
        else
            cost += T[tour[position-1]+1, tour[position+2]+1] - T[tour[position-1]+1, tour[position]+1] - T[tour[position]+1, tour[position+1]+1] - T[tour[position+1]+1, tour[position+2]+1]
        end
    end
    return cost
end

function calculate_new_cost_remove_three(tour::Vector{Int}, cost::Float64, position::Int, T::Matrix{Float64}, n_nodes::Int)
    nt = length(tour)
    if nt == 3
        return 0.0
    else
        if position == 1
            cost += T[1, tour[4]+1] - T[tour[1]+1, tour[2]+1] - T[tour[2]+1, tour[3]+1] - T[tour[3]+1, tour[4]+1] - T[1, tour[1]+1]
        elseif position == nt - 2
            cost += T[tour[nt-3]+1, n_nodes+2] - T[tour[nt-3]+1, tour[nt-2]+1] - T[tour[nt-2]+1, tour[nt-1]+1] - T[tour[nt-1]+1, tour[nt]+1] - T[tour[nt]+1, n_nodes+2]
        else
            cost += T[tour[position-1]+1, tour[position+3]+1] - T[tour[position-1]+1, tour[position]+1] - T[tour[position]+1, tour[position+1]+1] - T[tour[position+1]+1, tour[position+2]+1] - T[tour[position+2]+1, tour[position+3]+1]
        end
    end
    return cost
end

function calculate_new_cost_swap_one(tour1::Vector{Int}, depot1::Int64, cost1::Float64, city1::Int, position1::Int, tour2::Vector{Int}, depot2::Int64, cost2::Float64, city2::Int, position2::Int, T::Matrix{Float64})
    nt1 = length(tour1)
    nt2 = length(tour2)
    new_cost1 = cost1
    new_cost2 = cost2

    if nt1 == 1
        new_cost1 = T[depot1, city2]
    else
        if position1 == 1
            new_cost1 += T[depot1, city2] + T[city2, tour1[2]] - T[depot1, city1] - T[city1, tour1[2]]
        elseif position1 == nt1
            new_cost1 += T[tour1[nt1-1], city2] - T[tour1[nt1-1], city1]
        else
            new_cost1 += T[tour1[position1-1], city2] + T[city2, tour1[position1+1]] - T[tour1[position1-1], city1] - T[city1, tour1[position1+1]]
        end
    end

    if nt2 == 1
        new_cost2 = T[depot2, city1]
    else
        if position2 == 1
            new_cost2 += T[depot2, city1] + T[city1, tour2[2]] - T[depot2, city2] - T[city2, tour2[2]]
        elseif position2 == nt2
            new_cost2 += T[tour2[nt2-1], city1] - T[tour2[nt2-1], city2]
        else
            new_cost2 += T[tour2[position2-1], city1] + T[city1, tour2[position2+1]] - T[tour2[position2-1], city2] - T[city2, tour2[position2+1]]
        end
    end
    return new_cost1, new_cost2
end


function calculate_new_cost_swap_two_updated(tour1::Vector{Int}, cost1::Float64, city11::Int, city12::Int, position1::Int, tour2::Vector{Int}, cost2::Float64, city21::Int, city22::Int, position2::Int, T::Matrix{Float64}, n_nodes::Int)
    nt1 = length(tour1)
    nt2 = length(tour2)
    new_cost11 = cost1
    new_cost12 = cost1
    new_cost21 = cost2
    new_cost22 = cost2
    if nt1 == 2
        new_cost11 = T[1, city21+1] + T[city21+1, city22+1] + T[city22+1, n_nodes+2]
        new_cost12 = T[1, city22+1] + T[city22+1, city21+1] + T[city21+1, n_nodes+2]
    else
        if position1 == 1
            new_cost11 += T[1, city21+1] + T[city21+1, city22+1] + T[city22+1, tour1[3]+1] - T[1, city11+1] - T[city11+1, city12+1] - T[city12+1, tour1[3]+1]
            new_cost12 += T[1, city22+1] + T[city22+1, city21+1] + T[city21+1, tour1[3]+1] - T[1, city11+1] - T[city11+1, city12+1] - T[city12+1, tour1[3]+1]
        elseif position1 == nt1 - 1
            new_cost11 += T[tour1[nt1-2]+1, city21+1] + T[city21+1, city22+1] + T[city22+1, n_nodes+2] - T[tour1[nt1-2]+1, city11+1] - T[city11+1, city12+1] - T[city12+1, n_nodes+2]
            new_cost12 += T[tour1[nt1-2]+1, city22+1] + T[city22+1, city21+1] + T[city21+1, n_nodes+2] - T[tour1[nt1-2]+1, city11+1] - T[city11+1, city12+1] - T[city12+1, n_nodes+2]
        else
            new_cost11 += T[tour1[position1-1]+1, city21+1] + T[city21+1, city22+1] + T[city22+1, tour1[position1+2]+1] - T[tour1[position1-1]+1, city11+1] - T[city11+1, city12+1] - T[city12+1, tour1[position1+2]+1]
            new_cost12 += T[tour1[position1-1]+1, city22+1] + T[city22+1, city21+1] + T[city21+1, tour1[position1+2]+1] - T[tour1[position1-1]+1, city11+1] - T[city11+1, city12+1] - T[city12+1, tour1[position1+2]+1]
        end
    end

    if nt2 == 2
        new_cost21 = T[1, city11+1] + T[city11+1, city12+1] + T[city12+1, n_nodes+2]
        new_cost22 = T[1, city12+1] + T[city12+1, city11+1] + T[city11+1, n_nodes+2]
    else
        if position2 == 1
            new_cost21 += T[1, city11+1] + T[city11+1, city12+1] + T[city12+1, tour2[3]+1] - T[1, city21+1] - T[city21+1, city22+1] - T[city22+1, tour2[3]+1]
            new_cost22 += T[1, city12+1] + T[city12+1, city11+1] + T[city11+1, tour2[3]+1] - T[1, city21+1] - T[city21+1, city22+1] - T[city22+1, tour2[3]+1]
        elseif position2 == nt2 - 1
            new_cost21 += T[tour2[nt2-2]+1, city11+1] + T[city11+1, city12+1] + T[city12+1, n_nodes+2] - T[tour2[nt2-2]+1, city21+1] - T[city21+1, city22+1] - T[city22+1, n_nodes+2]
            new_cost22 += T[tour2[nt2-2]+1, city12+1] + T[city12+1, city11+1] + T[city11+1, n_nodes+2] - T[tour2[nt2-2]+1, city21+1] - T[city21+1, city22+1] - T[city22+1, n_nodes+2]
        else
            new_cost21 += T[tour2[position2-1]+1, city11+1] + T[city11+1, city12+1] + T[city12+1, tour2[position2+2]+1] - T[tour2[position2-1]+1, city21+1] - T[city21+1, city22+1] - T[city22+1, tour2[position2+2]+1]
            new_cost22 += T[tour2[position2-1]+1, city12+1] + T[city12+1, city11+1] + T[city11+1, tour2[position2+2]+1] - T[tour2[position2-1]+1, city21+1] - T[city21+1, city22+1] - T[city22+1, tour2[position2+2]+1]
        end
    end
    if new_cost11 > new_cost12
        if new_cost21 > new_cost22
            return new_cost11, new_cost21, true, true
        else
            return new_cost11, new_cost22, true, false
        end
    else
        if new_cost21 > new_cost22
            return new_cost12, new_cost21, false, true
        else
            return new_cost12, new_cost22, false, false
        end
    end
end

function calculate_new_cost_swap_three_updated(tour1::Vector{Int}, cost1::Float64, city11::Int, city12::Int, city13::Int, position1::Int, tour2::Vector{Int}, cost2::Float64, city21::Int, city22::Int, city23::Int, position2::Int, T::Matrix{Float64}, n_nodes::Int)
    nt1 = length(tour1)
    nt2 = length(tour2)
    new_cost11 = cost1
    new_cost12 = cost1
    new_cost21 = cost2
    new_cost22 = cost2
    if nt1 == 3
        new_cost11 = T[1, city21+1] + T[city21+1, city22+1] + T[city22+1, city23+1] + T[city23+1, n_nodes+2]
        new_cost12 = T[1, city23+1] + T[city23+1, city22+1] + T[city22+1, city21+1] + T[city21+1, n_nodes+2]
    else
        if position1 == 1
            new_cost11 += T[1, city21+1] + T[city21+1, city22+1] + T[city22+1, city23+1] + T[city23+1, tour1[4]+1] - T[1, city11+1] - T[city11+1, city12+1] - T[city12+1, city13+1] - T[city13+1, tour1[4]+1]
            new_cost12 += T[1, city23+1] + T[city23+1, city22+1] + T[city22+1, city21+1] + T[city21+1, tour1[4]+1] - T[1, city11+1] - T[city11+1, city12+1] - T[city12+1, city13+1] - T[city13+1, tour1[4]+1]
        elseif position1 == nt1 - 2
            new_cost11 += T[tour1[nt1-3]+1, city21+1] + T[city21+1, city22+1] + T[city22+1, city23+1] + T[city23+1, n_nodes+2] - T[tour1[nt1-3]+1, city11+1] - T[city11+1, city12+1] - T[city12+1, city13+1] - T[city13+1, n_nodes+2]
            new_cost12 += T[tour1[nt1-3]+1, city23+1] + T[city23+1, city22+1] + T[city22+1, city21+1] + T[city21+1, n_nodes+2] - T[tour1[nt1-3]+1, city11+1] - T[city11+1, city12+1] - T[city12+1, city13+1] - T[city13+1, n_nodes+2]
        else
            new_cost11 += T[tour1[position1-1]+1, city21+1] + T[city21+1, city22+1] + T[city22+1, city23+1] + T[city23+1, tour1[position1+3]+1] - T[tour1[position1-1]+1, city11+1] - T[city11+1, city12+1] - T[city12+1, city13+1] - T[city13+1, tour1[position1+3]+1]
            new_cost12 += T[tour1[position1-1]+1, city23+1] + T[city23+1, city22+1] + T[city22+1, city21+1] + T[city21+1, tour1[position1+3]+1] - T[tour1[position1-1]+1, city11+1] - T[city11+1, city12+1] - T[city12+1, city13+1] - T[city13+1, tour1[position1+3]+1]
        end
    end

    if nt2 == 3
        new_cost21 = T[1, city11+1] + T[city11+1, city12+1] + T[city12+1, city13+1] + T[city13+1, n_nodes+2]
        new_cost22 = T[1, city13+1] + T[city13+1, city12+1] + T[city12+1, city11+1] + T[city11+1, n_nodes+2]
    else
        if position2 == 1
            new_cost21 += T[1, city11+1] + T[city11+1, city12+1] + T[city12+1, city13+1] + T[city13+1, tour2[4]+1] - T[1, city21+1] - T[city21+1, city22+1] - T[city22+1, city23+1] - T[city23+1, tour2[4]+1]
            new_cost22 += T[1, city13+1] + T[city13+1, city12+1] + T[city12+1, city11+1] + T[city11+1, tour2[4]+1] - T[1, city21+1] - T[city21+1, city22+1] - T[city22+1, city23+1] - T[city23+1, tour2[4]+1]
        elseif position2 == nt2 - 2
            new_cost21 += T[tour2[nt2-3]+1, city11+1] + T[city11+1, city12+1] + T[city12+1, city13+1] + T[city13+1, n_nodes+2] - T[tour2[nt2-3]+1, city21+1] - T[city21+1, city22+1] - T[city22+1, city23+1] - T[city23+1, n_nodes+2]
            new_cost22 += T[tour2[nt2-3]+1, city13+1] + T[city13+1, city12+1] + T[city12+1, city11+1] + T[city11+1, n_nodes+2] - T[tour2[nt2-3]+1, city21+1] - T[city21+1, city22+1] - T[city22+1, city23+1] - T[city23+1, n_nodes+2]
        else
            new_cost21 += T[tour2[position2-1]+1, city11+1] + T[city11+1, city12+1] + T[city12+1, city13+1] + T[city13+1, tour2[position2+3]+1] - T[tour2[position2-1]+1, city21+1] - T[city21+1, city22+1] - T[city22+1, city23+1] - T[city23+1, tour2[position2+3]+1]
            new_cost22 += T[tour2[position2-1]+1, city13+1] + T[city13+1, city12+1] + T[city12+1, city11+1] + T[city11+1, tour2[position2+3]+1] - T[tour2[position2-1]+1, city21+1] - T[city21+1, city22+1] - T[city22+1, city23+1] - T[city23+1, tour2[position2+3]+1]
        end
    end
    if new_cost11 > new_cost12
        if new_cost21 > new_cost22
            return new_cost11, new_cost21, true, true
        else
            return new_cost11, new_cost22, true, false
        end
    else
        if new_cost21 > new_cost22
            return new_cost12, new_cost21, false, true
        else
            return new_cost12, new_cost22, false, false
        end
    end
end


function calculate_new_cost_swap_three_with_two_updated(tour1::Vector{Int}, cost1::Float64, city11::Int, city12::Int, city13::Int, position1::Int, tour2::Vector{Int}, cost2::Float64, city21::Int, city22::Int, position2::Int, T::Matrix{Float64}, n_nodes::Int)
    nt1 = length(tour1)
    nt2 = length(tour2)
    new_cost11 = cost1
    new_cost12 = cost1
    new_cost21 = cost2
    new_cost22 = cost2
    if nt1 == 3
        new_cost11 = T[1, city21+1] + T[city21+1, city22+1] + T[city22+1, n_nodes+2]
        new_cost12 = T[1, city22+1] + T[city22+1, city21+1] + T[city21+1, n_nodes+2]
    else
        if position1 == 1
            new_cost11 += T[1, city21+1] + T[city21+1, city22+1] + T[city22+1, tour1[4]+1] - T[1, city11+1] - T[city11+1, city12+1] - T[city12+1, city13+1] - T[city13+1, tour1[4]+1]
            new_cost12 += T[1, city22+1] + T[city22+1, city21+1] + T[city21+1, tour1[4]+1] - T[1, city11+1] - T[city11+1, city12+1] - T[city12+1, city13+1] - T[city13+1, tour1[4]+1]
        elseif position1 == nt1 - 2
            new_cost12 += T[tour1[nt1-3]+1, city22+1] + T[city22+1, city21+1] + T[city21+1, n_nodes+2] - T[tour1[nt1-3]+1, city11+1] - T[city11+1, city12+1] - T[city12+1, city13+1] - T[city13+1, n_nodes+2]
            new_cost11 += T[tour1[nt1-3]+1, city21+1] + T[city21+1, city22+1] + T[city22+1, n_nodes+2] - T[tour1[nt1-3]+1, city11+1] - T[city11+1, city12+1] - T[city12+1, city13+1] - T[city13+1, n_nodes+2]
        else
            new_cost11 += T[tour1[position1-1]+1, city21+1] + T[city21+1, city22+1] + T[city22+1, tour1[position1+3]+1] - T[tour1[position1-1]+1, city11+1] - T[city11+1, city12+1] - T[city12+1, city13+1] - T[city13+1, tour1[position1+3]+1]
            new_cost12 += T[tour1[position1-1]+1, city22+1] + T[city22+1, city21+1] + T[city21+1, tour1[position1+3]+1] - T[tour1[position1-1]+1, city11+1] - T[city11+1, city12+1] - T[city12+1, city13+1] - T[city13+1, tour1[position1+3]+1]
        end
    end

    if nt2 == 2
        new_cost21 = T[1, city11+1] + T[city11+1, city12+1] + T[city12+1, city13+1] + T[city13+1, n_nodes+2]
        new_cost22 = T[1, city13+1] + T[city13+1, city12+1] + T[city12+1, city11+1] + T[city11+1, n_nodes+2]
    else
        if position2 == 1
            new_cost21 += T[1, city11+1] + T[city11+1, city12+1] + T[city12+1, city13+1] + T[city13+1, tour2[3]+1] - T[1, city21+1] - T[city21+1, city22+1] - T[city22+1, tour2[3]+1]
            new_cost22 += T[1, city13+1] + T[city13+1, city12+1] + T[city12+1, city11+1] + T[city11+1, tour2[3]+1] - T[1, city21+1] - T[city21+1, city22+1] - T[city22+1, tour2[3]+1]
        elseif position2 == nt2 - 1
            new_cost21 += T[tour2[nt2-2]+1, city11+1] + T[city11+1, city12+1] + T[city12+1, city13+1] + T[city13+1, n_nodes+2] - T[tour2[nt2-2]+1, city21+1] - T[city21+1, city22+1] - T[city22+1, n_nodes+2]
            new_cost22 += T[tour2[nt2-2]+1, city13+1] + T[city13+1, city12+1] + T[city12+1, city11+1] + T[city11+1, n_nodes+2] - T[tour2[nt2-2]+1, city21+1] - T[city21+1, city22+1] - T[city22+1, n_nodes+2]
        else
            new_cost21 += T[tour2[position2-1]+1, city11+1] + T[city11+1, city12+1] + T[city12+1, city13+1] + T[city13+1, tour2[position2+2]+1] - T[tour2[position2-1]+1, city21+1] - T[city21+1, city22+1] - T[city22+1, tour2[position2+2]+1]
            new_cost22 += T[tour2[position2-1]+1, city13+1] + T[city13+1, city12+1] + T[city12+1, city11+1] + T[city11+1, tour2[position2+2]+1] - T[tour2[position2-1]+1, city21+1] - T[city21+1, city22+1] - T[city22+1, tour2[position2+2]+1]
        end
    end
    if new_cost11 > new_cost12
        if new_cost21 > new_cost22
            return new_cost11, new_cost21, true, true
        else
            return new_cost11, new_cost22, true, false
        end
    else
        if new_cost21 > new_cost22
            return new_cost12, new_cost21, false, true
        else
            return new_cost12, new_cost22, false, false
        end
    end
end


function calculate_new_cost_exchange_one(tour::Vector{Int}, depot::Int64, cost::Float64, city::Int, position1::Int, position2::Int, T::Matrix{Float64}, n_nodes::Int)
    nt = length(tour)
    if position1 == position2
        return cost
    end
    if position1 == 1
        cost = cost - T[depot, city] - T[city, tour[2]] + T[depot, tour[2]]
    elseif position1 == nt
        cost = cost - T[tour[nt-1], city]
    else
        cost = cost - T[tour[position1-1], city] - T[city, tour[position1+1]] + T[tour[position1-1], tour[position1+1]]
    end

    if position2 > position1
        if position2 == nt
            cost = cost + T[tour[nt], city]
        else
            cost = cost + T[tour[position2], city] + T[city, tour[position2+1]] - T[tour[position2], tour[position2+1]]
        end
    else
        if position2 == 1
            cost = cost + T[depot, city] + T[tour[1], city] - T[depot, tour[1]]
        else
            cost = cost + T[tour[position2-1], city] + T[city, tour[position2]] - T[tour[position2-1], tour[position2]]
        end
    end
    return cost
end

function calculate_new_cost_cross(tour1::Vector{Int}, cost1::Float64, tour2::Vector{Int}, cost2::Float64, k11::Int, k12::Int, k21::Int, k22::Int, T::Matrix{Float64}, n_nodes::Int)
    c1 = sum(T[tour1[i]+1, tour1[i+1]+1] for i in k11:k12-1)
    c2 = sum(T[tour2[i]+1, tour2[i+1]+1] for i in k21:k22-1)
    t1 = copy(tour1)
    t2 = copy(tour2)
    pushfirst!(t1, 0)
    push!(t1, n_nodes + 1)
    pushfirst!(t2, 0)
    push!(t2, n_nodes + 1)
    k11 += 1
    k12 += 1
    k21 += 1
    k22 += 1
    cost11 = cost1 - T[t1[k11-1]+1, t1[k11]+1] - c1 - T[t1[k12]+1, t1[k12+1]+1] + T[t1[k11-1]+1, t2[k21]+1] + c2 + T[t2[k22]+1, t1[k12+1]+1]
    cost12 = cost1 - T[t1[k11-1]+1, t1[k11]+1] - c1 - T[t1[k12]+1, t1[k12+1]+1] + T[t1[k11-1]+1, t2[k22]+1] + c2 + T[t2[k21]+1, t1[k12+1]+1]
    cost21 = cost2 - T[t2[k21-1]+1, t2[k21]+1] - c2 - T[t2[k22]+1, t2[k22+1]+1] + T[t2[k21-1]+1, t1[k11]+1] + c1 + T[t1[k12]+1, t2[k22+1]+1]
    cost22 = cost2 - T[t2[k21-1]+1, t2[k21]+1] - c2 - T[t2[k22]+1, t2[k22+1]+1] + T[t2[k21-1]+1, t1[k12]+1] + c1 + T[t1[k11]+1, t2[k22+1]+1]
    if cost11 < cost12
        if cost21 < cost22
            return cost11, cost21, true, true
        else
            return cost11, cost22, true, false
        end
    else
        if cost21 < cost22
            return cost12, cost21, false, true
        else
            return cost12, cost22, false, false
        end
    end
end

function calculate_new_cost_cross_upgraded(tour1::Vector{Int}, cost1::Float64, tour2::Vector{Int}, cost2::Float64, k11::Int, k12::Int, k21::Int, k22::Int, k13::Int, k23::Int, T::Matrix{Float64}, n_nodes::Int)
    c1 = sum(T[tour1[i]+1, tour1[i+1]+1] for i in k11:k12-1)
    c2 = sum(T[tour2[i]+1, tour2[i+1]+1] for i in k21:k22-1)
    t1 = copy(tour1)
    t2 = copy(tour2)
    pushfirst!(t1, 0)
    push!(t1, n_nodes + 1)
    pushfirst!(t2, 0)
    push!(t2, n_nodes + 1)
    k11 += 1
    k12 += 1
    k21 += 1
    k22 += 1
    k13 += 1
    k23 += 1
    if k13 == k11
        cost11 = cost1 - T[t1[k11-1]+1, t1[k11]+1] - c1 - T[t1[k12]+1, t1[k12+1]+1] + T[t1[k11-1]+1, t2[k21]+1] + c2 + T[t2[k22]+1, t1[k12+1]+1]
        cost12 = cost1 - T[t1[k11-1]+1, t1[k11]+1] - c1 - T[t1[k12]+1, t1[k12+1]+1] + T[t1[k11-1]+1, t2[k22]+1] + c2 + T[t2[k21]+1, t1[k12+1]+1]
    else
        cost11 = cost1 - T[t1[k11-1]+1, t1[k11]+1] - c1 - T[t1[k12]+1, t1[k12+1]+1] + T[t1[k11-1]+1, t1[k12+1]+1] + T[t1[k13-1]+1, t2[k21]+1] + c2 + T[t2[k22]+1, t1[k13]+1] - T[t1[k13-1]+1, t1[k13]+1]
        cost12 = cost1 - T[t1[k11-1]+1, t1[k11]+1] - c1 - T[t1[k12]+1, t1[k12+1]+1] + T[t1[k11-1]+1, t1[k12+1]+1] + T[t1[k13-1]+1, t2[k22]+1] + c2 + T[t2[k21]+1, t1[k13]+1] - T[t1[k13-1]+1, t1[k13]+1]
    end

    if k23 == k21
        cost21 = cost2 - T[t2[k21-1]+1, t2[k21]+1] - c2 - T[t2[k22]+1, t2[k22+1]+1] + T[t2[k21-1]+1, t1[k11]+1] + c1 + T[t1[k12]+1, t2[k22+1]+1]
        cost22 = cost2 - T[t2[k21-1]+1, t2[k21]+1] - c2 - T[t2[k22]+1, t2[k22+1]+1] + T[t2[k21-1]+1, t1[k12]+1] + c1 + T[t1[k11]+1, t2[k22+1]+1]
    else
        cost21 = cost2 - T[t2[k21-1]+1, t2[k21]+1] - c2 - T[t2[k22]+1, t2[k22+1]+1] + T[t2[k21-1]+1, t2[k22+1]+1] + T[t2[k23-1]+1, t1[k11]+1] + c1 + T[t1[k12]+1, t2[k23]+1] - T[t2[k23-1]+1, t2[k23]+1]
        cost22 = cost2 - T[t2[k21-1]+1, t2[k21]+1] - c2 - T[t2[k22]+1, t2[k22+1]+1] + T[t2[k21-1]+1, t2[k22+1]+1] + T[t2[k23-1]+1, t1[k12]+1] + c1 + T[t1[k11]+1, t2[k23]+1] - T[t2[k23-1]+1, t2[k23]+1]
    end
    if cost11 < cost12
        if cost21 < cost22
            return cost11, cost21, true, true
        else
            return cost11, cost22, true, false
        end
    else
        if cost21 < cost22
            return cost12, cost21, false, true
        else
            return cost12, cost22, false, false
        end
    end
end


function calculate_new_cost_exchange_two(tour::Vector{Int}, depot::Int64, cost::Float64, city1::Int, position1::Int, city2::Int,
    position2::Int, T::Matrix{Float64}, service::Vector{Int})
    nt = length(tour)

    if position1 == position2
        return cost
    end
    if nt == 2
        cost = T[depot, tour[2]] + service[tour[2]] + T[tour[2], tour[1]] + service[tour[1]]
        return cost
    end
    if position2 == position1 + 1
        if position1 == 1
            cost = cost - T[depot, city1] - T[city2, tour[3]] + T[depot, city2] + T[city1, tour[3]]
        elseif position2 == nt
            cost = cost - T[tour[nt-2], city1] + T[tour[nt-2], city2]
        else
            cost = cost - T[tour[position1-1], city1] - T[city2, tour[position2+1]] + T[tour[position1-1], city2] + T[city1, tour[position2+1]]
        end
        return cost
    end
    if position1 == position2 + 1
        if position2 == 1
            cost = cost - T[depot, city2] - T[city1, tour[3]] + T[depot, city1] + T[city2, tour[3]]
        elseif position1 == nt
            cost = cost - T[tour[nt-2], city2] + T[tour[nt-2], city1]
        else
            cost = cost - T[tour[position2-1], city2] - T[city1, tour[position1+1]] + T[tour[position2-1], city1] + T[city2, tour[position1+1]]
        end
        return cost
    end
    if position1 == 1
        cost = cost - T[depot, city1] - T[city1, tour[2]] + T[depot, tour[2]]
        cost = cost + T[depot, city2] + T[city2, tour[2]] - T[depot, tour[2]]
    elseif position1 == nt
        cost = cost - T[tour[nt-1], city1]
        cost = cost + T[tour[nt-1], city2]
    else
        cost = cost + T[tour[position1-1], city2] + T[city2, tour[position1+1]] - T[tour[position1-1], tour[position1+1]]
        cost = cost - T[tour[position1-1], city1] - T[city1, tour[position1+1]] + T[tour[position1-1], tour[position1+1]]
    end

    if position2 == 1
        cost = cost - T[depot, city2] - T[city2, tour[2]] + T[depot, tour[2]]
        cost = cost + T[depot, city1] + T[city1, tour[2]] - T[depot, tour[2]]
    elseif position2 == nt
        cost = cost - T[tour[nt-1], city2]
        cost = cost + T[tour[nt-1], city1]
    else
        cost = cost + T[tour[position2-1], city1] + T[city1, tour[position2+1]] - T[tour[position2-1], tour[position2+1]]
        cost = cost - T[tour[position2-1], city2] - T[city2, tour[position2+1]] + T[tour[position2-1], tour[position2+1]]
    end
    return cost
end

function calculate_new_cost_or_opt2(tour::Vector{Int}, depot::Int64, cost::Float64, city1::Int, position1::Int,
    city2::Int, position2::Int, T::Matrix{Float64}, n_nodes::Int)

    if position1 == position2
        return cost
    end

    nt = length(tour)
    if position1 == 1
        cost = cost - T[depot, city1] - T[city2, tour[3]] + T[depot, tour[3]]
    elseif position1 == nt - 1
        cost = cost - T[tour[nt-2], city1]
    else
        cost = cost - T[tour[position1-1], city1] - T[city2, tour[position1+2]] + T[tour[position1-1], tour[position1+2]]
    end

    if position1 < position2
        if position2 == nt - 1
            cost = cost + T[tour[nt], city1]
        else 
            cost = cost + T[tour[position2+1], city1] + T[city2, tour[position2+2]] - T[tour[position2+1], tour[position2+2]]
        end
    else
        if position2 == 1
            cost = cost + T[depot, city1] + T[city2, tour[1]] - T[depot, tour[1]]
        else
            cost = cost + T[tour[position2-1], city1] + T[city2, tour[position2]] - T[tour[position2-1], tour[position2]]
        end
    end
    return cost
end


function calculate_new_cost_or_opt3(tour::Vector{Int}, depot::Int64, cost::Float64, city1::Int, city2::Int, city3::Int, position1::Int,
    position2::Int, T::Matrix{Float64}, n_nodes::Int)

    if position1 == position2
        return cost
    end
    nt = length(tour)
    if position1 == 1
        cost = cost - T[depot, city1] - T[city3, tour[4]] + T[depot, tour[4]]
    elseif position1 == nt - 2
        cost = cost - T[tour[nt-3], city1]
    else
        cost = cost - T[tour[position1-1], city1] - T[city3, tour[position1+3]] + T[tour[position1-1], tour[position1+3]]
    end
    if position2 < position1
        if position2 == 1
            cost = cost + T[depot, city1] + T[city3, tour[1]] - T[depot, tour[1]]
        else
            cost = cost + T[tour[position2-1], city1] + T[city3, tour[position2]] - T[tour[position2-1], tour[position2]]
        end
    else
        if position2 == nt - 2
            cost = cost + T[tour[nt], city1]
        else
            cost = cost + T[tour[position2+2], city1] + T[city3, tour[position2+3]] - T[tour[position2+2], tour[position2+3]]
        end
    end
    return cost
end

function calculate_new_cost_2_opt(tour::Vector{Int}, depot::Int64, time::Float64, position1::Int, position2::Int, T::Matrix{Float64}, n_nodes::Int)

    nt = length(tour)
    if position1 == 1
        time = time - T[depot, tour[1]] + T[depot, tour[position2]]
    else
        time = time - T[tour[position1-1], tour[position1]] + T[tour[position1-1], tour[position2]]
    end

    if position2 == nt
        time = time
    else
        time = time - T[tour[position2], tour[position2+1]] + T[tour[position1], tour[position2+1]]
    end
    return time
end

function calculate_new_cost_3_opt(tour::Vector{Int}, depot::Int64, cost::Float64, k1::Int, k2::Int, k3::Int,
    T::Matrix{Float64}, n_nodes::Int)

    nt = length(tour)
    if k3 - k2 > 2
        position1 = k2 + 1
        position2 = k3 - 1
        cost = cost - T[tour[position1-1], tour[position1]] + T[tour[position1-1], tour[position2]] - T[tour[position2], tour[position2+1]] + T[tour[position1], tour[position2+1]]
    end
    if k2 - k1 > 2
        position1 = k1 + 1
        position2 = k2 - 1
        cost = cost - T[tour[position1-1], tour[position1]] + T[tour[position1-1], tour[position2]] - T[tour[position2], tour[position2+1]] + T[tour[position1], tour[position2+1]]
    end
    return cost
end

function calculate_new_cost_3_permute(tour::Vector{Int}, cost::Float64, S1::Vector{Int}, S2::Vector{Int}, k1::Int, T::Matrix{Float64}, n_nodes::Int)
    if S1 == S2
        return cost
    end
    nt = length(tour)
    if nt == 3
        cost = T[1, S2[1]+1] + T[S2[1]+1, S2[2]+1] + T[S2[2]+1, S2[3]+1] + T[S2[3]+1, n_nodes+2]
        return cost
    end
    cost = cost - T[S1[1]+1, S1[2]+1] - T[S1[2]+1, S1[3]+1] + T[S2[1]+1, S2[2]+1] + T[S2[2]+1, S2[3]+1]

    if k1 == 1
        cost = cost - T[1, S1[1]+1] + T[1, S2[1]+1] - T[S1[3]+1, tour[4]+1] + T[S2[3]+1, tour[4]+1]
    elseif k1 == nt - 2
        cost = cost - T[S1[3]+1, n_nodes+2] + T[S2[3]+1, n_nodes+2] - T[tour[nt-3]+1, S1[1]+1] + T[tour[nt-3]+1, S2[1]+1]
    else
        cost = cost - T[tour[k1-1]+1, S1[1]+1] + T[tour[k1-1]+1, S2[1]+1] - T[S1[3]+1, tour[k1+3]+1] + T[S2[3]+1, tour[k1+3]+1]
    end
    return cost
end