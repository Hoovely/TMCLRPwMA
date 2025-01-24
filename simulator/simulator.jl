mutable struct env
    routes::Vector{Int64}
    weighted_time::Float64
    g_obj::Float64
    power::Float64
    total_cost::Float64
    depots::Vector{Int64}
    tours::Vector{Tour}
    reds_path::Vector{Tour}
    r_cost::Float64
    r_obj::Float64
end

include("greedy_policy.jl")
include("set_covering_policy.jl")



function divided_budget_policy(TMC_num::Int64, green_num::Int64, red_num::Int64, EC_G::Int64, EC_R::Int64, VC_G::Int64, VC_R::Int64, FC::Int64, epsilon::Int64)
    
    exp_cost = 0
    exp_cost += red_num * (VC_R+EC_R) * 1.1
    
    number_g_ambulance = floor(Int64, green_num/2)
    exp_cost += number_g_ambulance * (VC_G+EC_G) * 0.5

    construction_cost = epsilon - exp_cost
    opening_facility = floor(Int64, construction_cost/FC)

    return opening_facility, number_g_ambulance
end

function simulator(problem::String, node_size::Int64, TMC_num::Int64, green_num::Int64, red_num::Int64, hospital_num::Int64, WG::Int64, WR::Int64, 
    service::Vector{Int64}, patient::Vector{Int64}, demand::Vector{Int64}, coordinate::Vector{Vector}, distance::Matrix{Float64}, travel::Matrix{Float64},
    N::Vector{Int64}, F::Vector{Int64}, G::Vector{Int64}, R::Vector{Int64}, H::Vector{Int64}, return_tmc::Vector{Int64}, capacity_vg::Int64, capacity_vr::Int64, 
    capacity_room::Vector{Int64}, capacity_man::Vector{Int64}, EC_G::Int64, EC_R::Int64, VC_G::Int64, VC_R::Int64, FC::Int64, epsilon::Int64, instance_str::String, 
    simulator_version::String, policy::String, verbose::Bool)

    println("===========================================================================================================")
    println("Start ", policy, " Policy Simulation")

    t1 = time()

    num_open, num_non_ambulance = divided_budget_policy(TMC_num, green_num, red_num, EC_G, EC_R, VC_G, VC_R, FC, epsilon)
    
    if policy == "Greedy"
        F_ = greedy_solve_location_problem(num_open, WG, WR, distance, F, G, R, FC)
        assigned_cluster_g, assigned_cluster_r = greedy_solve_assignment_problem(distance, F_, G, R) 
        ne_tours, e_tours, obj, cost = greedy_solve_routing_problem(node_size, F_, H, assigned_cluster_g, assigned_cluster_r, WG, WR, 
                                                                    service, patient, demand, distance, travel, capacity_vg, capacity_man,
                                                                    EC_G, EC_R, VC_G, VC_R, FC)
    end

    if policy == "Set Covering"
        radius = 1000.0
        println(num_open)
        F_ = set_covering_solve_location_problem(num_open, travel, F, G, R, radius)
        println(F_)
        assigned_cluster_g, assigned_cluster_r = greedy_solve_assignment_problem(distance, F_, G, R) 
        ne_tours, e_tours, obj, cost = greedy_solve_routing_problem(node_size, F_, H, assigned_cluster_g, assigned_cluster_r, WG, WR, service, patient, demand, distance, travel, 
                                                             capacity_vg, capacity_man, EC_G, EC_R, VC_G, VC_R, FC)
    end

    calculation_time = time() - t1

    println("===========================================================================================================")
    println("Simulation Obj: ", round(obj, digits=3))
    println("Operation Cost: ", round(cost, digits=3))
    println("Total Run time: ", round(calculation_time, digits=3))
    
end