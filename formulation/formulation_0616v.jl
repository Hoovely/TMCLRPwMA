include("../data_generator/data_generator.jl")
include("../data_generator/augerat_matrix_generation_cut.jl")
include("../data_generator/load_data_augerat.jl")
include("../data_generator/duplicate_data_augerat.jl")

using JuMP, Gurobi, CPLEX, CSV

function formulation(budget::Int64, cutting_num::Int64, N::Vector{Int64}, F::Vector{Int64}, G::Vector{Int64}, R::Vector{Int64}, H::Vector{Int64}, K::Vector{Int64}, Return_TMC::Vector{Int64}, Service::Vector{Int64}, Travel::Matrix{Float64},
    patient::Vector{Int64}, Demand::Vector{Int64}, Capacity_V::Int64, Capacity_H::Int64, EC_G::Int64, EC_R::Int64, VC_G::Int64, VC_R::Int64, FC::Int64, M::Int64, D::Int64, WG::Int64, WR::Int64)


    # Initialize Cplex model
    model = Model(CPLEX.Optimizer)

    global x, y, z, a, eG, eR, b, u
    x = @variable(model, x[i in N, j in N, k in K], Bin)
    y = @variable(model, y[f in F], Bin)
    z = @variable(model, z[i in vcat(G, R), f in F, k in K], Bin)
    a = @variable(model, a[f in F, k in K], Int)
    eG = @variable(model, eG, Int)
    eR = @variable(model, eR, Int)
    b = @variable(model, b[i in vcat(F, G, R)], Int)
    u = @variable(model, u[i in R, h in H], Bin)

    # Objective function
    @objective(model, Min, WR * eR + WG * eG)

    
    # Constraints
    # Epsilon Constraint
    # Constraint 2
    @constraint(model, EC_G * sum(a[f, k] for f in F for k in K) + 
                       VC_G * sum(x[f, g, k] for f in F for g in G for k in K) +
                       VC_R * sum(x[f, r, k] for f in F for r in R for k in K) + 
                       FC * sum(y[f] for f in F) <= budget )

    # Multi-Depot Constraints
    # Constraint 3
    for j in G
        @constraint(model, sum(x[i, j, k] for k in K for i in vcat(F, G) if i != j) == 1)
    end

    # Constraint 4
    for i in G
        @constraint(model, sum(x[i, j, k] for k in K for j in vcat(Return_TMC, G) if i != j) == 1)
    end

    # Constraint 5
    for k in K
        for j in G
            @constraint(model, sum(x[i, j, k] for i in vcat(F, G) if i != j) == sum(x[j, i, k] for i in vcat(Return_TMC, G) if i != j))
        end
    end

    # Constraint 6
    for k in K
        for f in F
            @constraint(model, sum(x[f, j, k] for j in G) == sum(x[i, f+cutting_num, k] for i in G))
        end
    end

    # Constraint 7
    for k in K
        @constraint(model, sum(x[f, j, k] for f in F for j in vcat(G, R)) <= 1)
    end

    # Constraint 8
    for k in K
        @constraint(model, sum(x[i, f, k] for f in Return_TMC for i in G) <= 1)
    end

    # Emergency Patient
    # Constraint 9
    for j in R
        @constraint(model, sum(x[f, j, k] for f in F for k in K) == 1)
    end

    # Constraint 10
    for k in K
        for i in R
            @constraint(model, sum(x[f, i, k] for f in F) == sum(x[i, h, k] for h in H))
        end
    end

    # Vehicle constraints
    # Constraint 11
    for i in vcat(G, R)
        @constraint(model, sum(z[i, f, k] for f in F for k in K) <= 1)
    end

    # Constraint 12
    for f in F
        @constraint(model, sum(z[i, f, k] for k in K for i in vcat(G, R)) <= Capacity_V)
    end

    # Facility Location
    # Constraint 13
    for f in F
        for k in K
            for i in vcat(G, R)
                @constraint(model, z[i, f, k] <= y[f])
            end
        end
    end

    # Constraint 14
    for f in F
        for k in K
            for i in vcat(G, R)
                @constraint(model, x[f, i, k] <= y[f])
            end
        end
    end

    # Constraint 15
    for f in F
        for k in K
            for i in G
                @constraint(model, x[i, f+cutting_num, k] <= y[f])
            end
        end
    end

    # Manpower Allocation
    # Constraint 16
    for f in F
        for k in K
            for i in vcat(G, R)
                @constraint(model, a[f, k] - z[i, f, k] * Demand[i] >= 0)
            end
        end
    end


    # Constraint 18
    for f in F
        for k in K
            for j in vcat(G, R)
                @constraint(model, sum(x[f, i, k] for i in vcat(G, R)) + sum(x[i, j, k] for i in N if i != j) <= 1 + z[j, f, k])
            end
        end
    end

    # Rescue Completion Time
    # Constraint 19
    for i in R
        @constraint(model, sum(u[i, h] for h in H) == 1)
    end

    # Constraint 20
    for h in H
        @constraint(model, sum(u[i, h] for i in R) <= Capacity_H)
    end

    # Constraint 21
    for i in R
        for h in H
            @constraint(model, sum(x[i, h, k] for k in K) <= u[i, h])
        end
    end

    # Constraint 22
    for i in F
        for j in vcat(G, R)
            @constraint(model, Travel[i, j] <= b[j] + D * (1 - sum(x[i, j, k] for k in K)))
        end
    end

    # Constraint 23
    for i in G
        for j in G
            if i != j
                @constraint(model, b[i] + Service[i] + Travel[i, j] <= b[j] + D * (1 - sum(x[i, j, k] for k in K)))
            end
        end
    end

    # Constraint 24
    for i in G
        @constraint(model, eG >= b[i] + Service[i])
    end

    # Constraint 25
    for i in R
        for h in H
            @constraint(model, eR >= b[i] + Service[i] + u[i, h] * Travel[i, h])
        end
    end

    return model
end


function solve_MIP(model, cpu, timelimit=3600)
    set_time_limit_sec(model, timelimit)
    set_optimizer_attribute(model, "CPXPARAM_Threads", cpu)

    start_time = time()
    optimize!(model)
    solve_time = time() - start_time

    return objective_value(model), MOI.get(model, MOI.RelativeGap()), solve_time

end


function callback_variable(N::Vector{Int64}, F::Vector{Int64}, G::Vector{Int64}, R::Vector{Int64}, H::Vector{Int64}, K::Vector{Int64})

    x_sol = [(i, j, k, value.(x[i, j, k])) for i in N for j in N for k in K if value.(x[i, j, k]) >= 0.9]
    y_sol = [(f, value.(y[f])) for f in F]
    z_sol = [(i, f, k, value.(z[i, f, k])) for i in vcat(G, R) for f in F for k in K if value.(z[i, f, k]) >= 0.9]
    a_sol = [(f, k, value.(a[f, k])) for f in F for k in K if value.(a[f, k]) >= 0.9]
    u_sol = [(i, h, value.(u[i, h])) for i in R for h in H if value.(u[i, h]) >= 0.9]
    b_sol = [(i, value.(b[i])) for i in vcat(F, G, R) if value.(b[i]) >= 0.9]
    eG_sol = value.(eG)
    eR_sol = value.(eR)
    
    return [x_sol, y_sol, z_sol, a_sol, u_sol, b_sol, eG_sol, eR_sol]
    
end

function draw_route(coord, problem, TMC_num, green_num, red_num, hospital_num, cutting_num, obj, epsilon, instance_str, log_dir)
    
    x_coordinate, y_coordinate = [], []
    for (x,y) in coord
        push!(x_coordinate, x)
        push!(y_coordinate, y)
    end
    
    # Plot
    plot()
    
    # Main Route
    for (i, j, k, val) in x_sol
        if val >= 0.9
            x1, y1 = x_coordinate[i], y_coordinate[i]
            x2, y2 = x_coordinate[j], y_coordinate[j]
            arrow!(plot(), x1, y1, x2, y2, line=:dash, color=:black)
        end
    end
    
    # Plot node
    for i in 1:length(x_coordinate)
        if i <= TMC_num
            if any(x -> x[1] == i && x[2] >= 0.9, y_sol)
                scatter!(plot(), [x_coordinate[i]], [y_coordinate[i]], color=:blue, label="Facility")
            else 
                scatter!(plot(), [x_coordinate[i]], [y_coordinate[i]], color=:black)
            end
        elseif i <= TMC_num + green_num
            scatter!(plot(), [x_coordinate[i]], [y_coordinate[i]], color=:green, label="Green Node")
        elseif i <= TMC_num + green_num + red_num
            scatter!(plot(), [x_coordinate[i]], [y_coordinate[i]], color=:red, label="Red Node")
        elseif i <= TMC_num + green_num + red_num + hospital_num
            scatter!(plot(), [x_coordinate[i]], [y_coordinate[i]], color=:purple, label="Hospital")
        end
    end

    title!("Optimal solution for $problem, instance setting: $instance_str, sol: $obj, $epsilon")
    log_dir = joinpath(log_dir, "figure")
    mkpath(log_dir)
    savefig(joinpath(log_dir, "$instance_str, $obj, $epsilon.png"))
end





function run_math(problem::String, coordinates::Vector{Vector{Int64}}, cutting_num::Int64, TMC_num::Int64, green_num::Int64, red_num::Int64, hospital_num::Int64, WG::Int64, WR::Int64, epsilon::Int64,
    N::Vector{Int64}, F::Vector{Int64}, G::Vector{Int64}, R::Vector{Int64}, H::Vector{Int64}, K::Vector{Int64}, return_tmc::Vector{Int64}, service::Vector{Int64}, travel::Matrix{Float64},
    patient::Vector{Int64}, demand::Vector{Int64}, capacity_v::Int64, capacity_h::Int64, EC_G::Int64, EC_R::Int64, VC_G::Int64, VC_R::Int64, FC::Int64, M::Int64, D::Int64, instance_str::String, formulation_version::String)
    

    model = formulation(epsilon, cutting_num, N, F, G, R, H, K, return_tmc, service, travel, patient, demand, capacity_v, capacity_h, EC_G, EC_R, VC_G, VC_R, FC, M, D, WG, WR)
    obj, gap, calculation_time = solve_MIP(model, 15)
    solution_list = callback_variable(N, F, G, R, H, K)

    df = DataFrame()
    df["Problem"] = [problem]
    df["Node Size"] = [cutting_num]
    df["Temporary Emergency Center Number"] = [TMC_num]
    df["Non-Emergency Patient Number"] = [green_num]
    df["Emergency Patient Number"] = [red_num]
    df["Hospital Number"] = [hospital_num]
    df["Non-Emergency Patient Weights"] = [WG]
    df["Emergency Patient Weights"] = [WR]
    df["Epsilon Value"] = [epsilon]
    df["Obj"] = [obj]
    df["CPU Time"] = [calculation_time]
    df["Relative Gap"] = [gap * 100]

    log_dir = "result/mathematical_programming"
    mkpath(log_dir)

    formulation_version = "v0616"
    log_dir = joinpath(log_dir, formulation_version)
    mkpath(log_dir)

    if !isdir(joinpath(log_dir, "solution"))
        mkpath(joinpath(log_dir, "solution"))
    end

    draw_route(coordinates, problem, TMC_num, green_num, red_num, hospital_num, cutting_num, obj, epsilon, instance_str, log_dir)

    csv_file = joinpath(log_dir, "TMCLRPwMA_fomulation_$formulation_version.csv")
    if !isfile(csv_file)
        CSV.write(csv_file, df)
    else
        CSV.write(csv_file, df, append=true)
    end

end


