using Clustering
using Distances
using StatsBase
using TSPLIB


# using Concorde
# using TSPSolvers
# using LKH
import TravelingSalesmanHeuristics: solve_tsp as TSP_Heuristic_solve_tsp
import Hygese: solve_tsp as Hygese_solve_tsp
using Hygese

# using LinearAlgebra


include("split.jl")
include("genetic_algorithm.jl")
include("initial.jl")
include("rebalancing.jl")
include("red_rescue.jl")
include("mutation.jl")
include("crossover.jl")
include("neighborhood.jl")
include("neighborhood_intra.jl")
include("costs.jl")
include("intersection.jl")
include("enrichment.jl")


function callback_variable(chrm::Chromosome)

    x_sol = Vector{Vector}()
    y_sol = Int[]

    for tour in chrm.tours
        depot = tour.sequence[1]
        push!(y_sol, depot)

        for (idx, node) in enumerate(tour.sequence)
            if node == tour.sequence[end]
                push!(x_sol, [node, depot])
                break
            end

            next_node = tour.sequence[idx+1]
            push!(x_sol, [node, next_node])
        end
    end

    red_paths = Vector{Vector}()
    for tour in chrm.reds_path
        push!(y_sol, tour.sequence[1])
        push!(x_sol, [tour.sequence[1], tour.sequence[2]])
        push!(x_sol, [tour.sequence[2], tour.sequence[3]])

        push!(red_paths, [tour.sequence[1], tour.sequence[2], tour.sequence[3]])
    end

    println(x_sol)

    y_sol = collect(Set(y_sol))

    return [x_sol, y_sol, red_paths]
    
end

function save_result(problem, cutting_num, TMC_num, green_num, red_num, hospital_num, WG, WR, epsilon, 
                     best_cost, best_obj, aver_cost, aver_obj, worst_cost, worst_obj, total_time, aver_time, heuristic_version, 
                     solution_list, coordinates, instance_str, case_study_flag)
    
    df = DataFrame(
        Problem=[problem],
        Node_Size=[cutting_num],
        Temporary_Emergency_Center_Number=[TMC_num],
        Non_Emergency_Patient_Number=[green_num],
        Emergency_Patient_Number = [red_num],
        Hospital_Number = [hospital_num],
        Non_Emergency_Patient_Weights = [WG],
        Emergency_Patient_Weights = [WR],
        Epsilon_Value = [epsilon],
        Best_Obj = [best_obj],
        Best_Cost = [best_cost],
        Average_Obj = [aver_obj],
        Average_Cost = [aver_cost],
        Worst_Obj = [worst_obj],
        Worst_Cost = [worst_cost],
        Time = [aver_time],
    )

    log_dir = "result/heuristic"
    mkpath(log_dir)

    log_dir = joinpath(log_dir, heuristic_version)
    mkpath(log_dir)

    if !isdir(joinpath(log_dir, "solution"))
        mkpath(joinpath(log_dir, "solution"))
    end

    csv_file = joinpath(log_dir, "Hybrid_GA_$heuristic_version.csv")
    if !isfile(csv_file)
        CSV.write(csv_file, df)
    else
        CSV.write(csv_file, df, append=true)
    end

    if case_study_flag
        red_paths = solution_list[3]

        x_coordinate, y_coordinate = [], []
        for (x,y) in coordinates
            push!(x_coordinate, x)
            push!(y_coordinate, y)
        end

        mkpath(joinpath(log_dir, "save_dispatch_paths"))
        file = open(string(joinpath(log_dir, "save_dispatch_paths"), instance_str), "w")

        for (f, r, h) in red_paths
            x1, y1 = x_coordinate[f], y_coordinate[f]
            x2, y2 = x_coordinate[r], y_coordinate[r]
            x3, y3 = x_coordinate[h], y_coordinate[h]
            
            println(file, "$(x1)\t$(y1)\t$(x2)\t$(y2)\t$(x3)\t$(y3)")
        end

        close(file)

    end

    draw_route(solution_list, coordinates, problem, TMC_num, green_num, red_num, hospital_num, cutting_num, best_obj, best_cost, epsilon, instance_str, log_dir) 

end

function draw_route(solution_list, coord, problem, TMC_num, green_num, red_num, hospital_num, cutting_num, obj, operation_cost, epsilon, instance_str, log_dir)
    x_sol, y_sol = solution_list[1], solution_list[2]

    x_coordinate, y_coordinate = [], []
    for (x,y) in coord
        push!(x_coordinate, x)
        push!(y_coordinate, y)
    end

    plot()

    # Main Route
    for (i, j) in x_sol
        x1, y1 = x_coordinate[i], y_coordinate[i]
        x2, y2 = x_coordinate[j], y_coordinate[j]
        
        plot!([x1, x2], [y1, y2], linewidth=1.0, color="black", legend=:none)
    end

    # Plot node
    for i in 1:length(x_coordinate)
        if i <= TMC_num
            if any(x -> x[1] == i, y_sol)
                scatter!((x_coordinate[i], y_coordinate[i]), color="blue", legend=:none)
            else 
                scatter!((x_coordinate[i], y_coordinate[i]), color="black", legend=:none)
            end
        elseif i <= TMC_num + green_num
            scatter!((x_coordinate[i], y_coordinate[i]), color="green", legend=:none)
        elseif i <= TMC_num + green_num + red_num
            scatter!((x_coordinate[i], y_coordinate[i]), color="red", legend=:none)
        elseif i <= TMC_num + green_num + red_num + hospital_num
            scatter!((x_coordinate[i], y_coordinate[i]), color="purple", legend=:none)
        end
    end

    title!("Solution of GA for $problem, instance setting: $instance_str, budget: $epsilon, sol: ($obj, $operation_cost)")
    log_dir = joinpath(log_dir, "figure")
    mkpath(log_dir)
    savefig(joinpath(log_dir, "$instance_str, $obj, $operation_cost.png"))

end

function draw_obj(obj_list, index, instance_str, heuristic_version)
    
    df = DataFrame(obj_list=obj_list)
    
    log_dir = "result/heuristic"
    mkpath(log_dir)

    log_dir = joinpath(log_dir, heuristic_version)
    mkpath(log_dir)

    if !isdir(joinpath(log_dir, "solution"))
        mkpath(joinpath(log_dir, "solution"))
    end

    log_dir = joinpath(log_dir, "obj_csv")
    mkpath(log_dir)
    CSV.write(joinpath(log_dir, "$index, $instance_str.csv"), df)

end

function run_GA(problem::String, node_size::Int64, TMC_num::Int64, green_num::Int64, red_num::Int64, hospital_num::Int64, WG::Int64, WR::Int64, 
    service::Vector{Int64}, patient::Vector{Int64}, demand::Vector{Int64}, coordinate::Vector{Vector}, distance::Matrix{Float64}, travel::Matrix{Float64},
    N::Vector{Int64}, F::Vector{Int64}, G::Vector{Int64}, R::Vector{Int64}, H::Vector{Int64}, return_tmc::Vector{Int64}, capacity_vg::Int64, capacity_vr::Int64, 
    capacity_room::Vector{Int64}, capacity_man::Vector{Int64}, EC_G::Int64, EC_R::Int64, VC_G::Int64, VC_R::Int64, FC::Int64, epsilon::Int64, instance_str::String, 
    heuristic_version::String, verbose::Bool, case_study_flag::Bool)

    ### Setting Problem Set, Parameter ###
    # F, G, R, H, Return_Depot = create_set(num_depot, num_green, num_red, num_hospital, N)
    # Patient, Demand, C, K, FC, EC = create_parameter(node_size, num_green, num_red, F, G, R)

    # K = trunc(Int64, ceil(Int64, sum(patient[1:length(F)+length(G)]) / (capacity_v-maximum(demand[1:length(F)+length(G)])) ))

    # randomly coordinates based
    # Customers = rand(length(G), 2) 

    # real coordinates based
    cus = coordinate[length(F)+1:length(F)+length(G)]
    Customers = Matrix{Float64}(undef, length(G), 2)
    for i in 1:length(G)
        Customers[i, 1] = cus[i][1]
        Customers[i, 2] = cus[i][2]
    end 

    depo = coordinate[1:length(F)]
    Depots = Matrix{Float64}(undef, length(F), 2)
    for i in 1:length(F)
        Depots[i, 1] = depo[i][1]
        Depots[i, 2] = depo[i][2]
    end

    # depot, Customers = create_random_sample(num_customers) # Randomly create coordinates sample
    # depot_coordinates, customer_coordinates = coordinates[1, :], coordinates[2:end, :] # Using dataset instance coordinates

    ### GA Parameter ###
    # n = size(T)[1] - 2
    h = 0.3
    popsize = (10, 20)
    k_tournament = 2
    num_iter = 250
    time_limit = 10000.0
    Mutation_Chance = 0.5
    num_runs = 10
    num_nei = 2
    crossover_functions = [2]

    
    ### Start GA ###
    P = Chromosome[]
    Elite_P = Chromosome[]

    avg_obj = 0.0
    best_obj = Inf
    worst_obj = 0.0

    avg_cost = 0.0
    best_cost = Inf
    worst_cost = 0.0

    t1 = time()
    for i in 1:num_runs
        println("===========================================================================================================")
        println("Iter $(i) Start!")

        P, roullet, obj_vec = perform_genetic_algorithm(node_size, travel, h, popsize, k_tournament, num_iter, time_limit, Mutation_Chance, num_nei, 
                                                crossover_functions, Customers, Depots, F, G, R, H, return_tmc, distance, service, patient, demand,
                                                capacity_vg, capacity_vr, capacity_room, capacity_man, FC, EC_G, EC_R, VC_G, VC_R, WG, WR, coordinate, epsilon, verbose)
        avg_obj += P.weighted_time
        avg_cost += P.total_cost
        if P.weighted_time < best_obj
            best_obj = P.weighted_time
            best_cost = P.total_cost
            Elite_P = P
        end
        if P.weighted_time > worst_obj
            worst_obj = P.weighted_time
            worst_cost = P.total_cost
        end

        println("Elapsed time: $( round(time()-t1, digits=3) ) sec")

        draw_obj(obj_vec, i, instance_str, heuristic_version)
    end
    calculation_time = time() - t1

    println("===========================================================================================================")
    println("Best Obj: ", round(best_obj, digits=3), "  Average Obj: ", round(avg_obj / num_runs, digits=3), "  Worst Obj: ", round(worst_obj, digits=3))
    println("Best Route Operation Cost: ", round(best_cost, digits=3), "  Average Route Operation Cost: ", round(avg_cost / num_runs, digits=3), "  Worst Route Operation Cost: ", round(worst_cost, digits=3))
    println("Total Run time: ", round(calculation_time, digits=3), ", Average Run time: ", round(calculation_time / num_runs, digits=3))
    best_route(Elite_P)
    
    solution_list = callback_variable(Elite_P)
    save_result(problem, node_size, TMC_num, green_num, red_num, hospital_num, WG, WR, epsilon, 
                best_cost, best_obj, avg_cost/num_runs, avg_obj/num_runs, worst_cost, worst_obj, calculation_time, calculation_time/num_runs, 
                heuristic_version, solution_list, coordinate, instance_str, case_study_flag)

end
