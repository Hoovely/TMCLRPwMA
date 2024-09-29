import Pandas.read_pickle

include("data_generator/prodhon_matrix_generation.jl")
include("formulation/formulation_0929v.jl")
include("heuristic/src/main.jl")


function load_instance(folder_path)
    files = filter(isfile, readdir(folder_path, join=true))

    data = Vector()
    for file in files
        instance = Vector()

        open(file, "r") do f
            for l in eachline(f)
                if l == "" continue end
                
                li = split(l, "\t")
                if length(li) == 1
                    lin = parse(Float64, li[1])
                    line = isinteger(lin) ? Int(lin) : lin
                else
                    line = Vector()
                    for c in li
                        if c == "" continue end
                        lin = parse(Float64, c)
                        push!(line, isinteger(lin) ? Int(lin) : lin)
                    end
                end

                push!(instance, line)
            end
        end
        push!(data, [file, instance])
    end

    instance_set = []
    for (file, problem) in data

        number_f = problem[1]
        number_g = problem[2]
        number_r = problem[3]
        number_h = problem[4]
    
        coordi = Vector{Vector}()
        i = 5
        j = 4+number_f+number_g+number_r+number_h
        for (x, y) in problem[i:j]
            push!(coordi, [x,y])
        end
    
        service = Vector{Int}()
        i = j+1
        j = j+number_f+number_g+number_r+number_h
        for v in problem[i:j]
            push!(service, v)
        end
    
        patient = Vector{Int}()
        i = j+1
        j = j+number_f+number_g+number_r+number_h
        for v in problem[i:j]
            push!(patient, v)
        end
    
        demand = Vector{Int}()
        i = j+1
        j = j+number_f+number_g+number_r+number_h
        for v in problem[i:j]
            push!(demand, v)
        end

        j = j+1
        budget = problem[j]

        capa_room = Vector{Int}()
        i = j+1
        j = j+number_h
        for v in problem[i:j]
            push!(capa_room, v)
        end

        capa_man = Vector{Int}()
        i = j+1
        j = j+number_h
        for v in problem[i:j]
            push!(capa_man, v)
        end
    
        push!(instance_set, [file, number_f, number_g, number_r, number_h, coordi, service, patient, demand, budget, capa_room, capa_man])
    end

    return instance_set
end


function generate_set_and_parameter(cutting_num::Int64, TMC_num::Int64, green_num::Int64, red_num::Int64, hospital_num::Int64)
    
    # Set parameters
    N = [i for i in 1:cutting_num+TMC_num]
    F = generate_node(N, TMC_num)
    G = generate_node(difference(N, F), green_num)
    R = generate_node(difference(N, vcat(F, G)), red_num)
    H = generate_node(difference(N, vcat(F, G, R)), hospital_num)
    K = collect(1:(green_num + red_num))
    Return_TMC = difference(N, vcat(F, G, R, H))

    # Capacities
    Capacity_VG, Capacity_VR = 12, 5

    # Costs
    EC_G = 30
    EC_R = 50
    VC_G = 1000
    VC_R = 1000
    FC = 20000

    # Big-M
    M = 1000000

    return N, F, G, R, H, K, Return_TMC, Capacity_VG, Capacity_VR, EC_G, EC_R, VC_G, VC_R, FC, M
end


function prodhon_experiment()
    folder_path = "data_generator/TMCLRP_Instance"
    instance_set = load_instance(folder_path)

    for (WG, WR) in [(1,1)]
        for (problem, number_f, number_g, number_r, number_h, coordi, service, patient, demand, budget, capa_room, capa_man) in instance_set    
            if number_g + number_r >= 100 continue end

            cutting_num = number_f + number_g + number_r + number_h
            instance_str = "($cutting_num, $number_f, $number_g, $number_r, $number_h, $budget)"
            distance_matrix = prodhon_matrix(coordi, number_f)
            Travel = generate_travel(distance_matrix, 1)
            set_and_parameter = generate_set_and_parameter(cutting_num, number_f, number_g, number_r, number_h)

            println("CURRENTLY SOLVING ==> $problem")
            println("NODE SIZE ==> $cutting_num")
            println("TEMPORARY MEDICAL CENTER NUMBER ==> $number_f")
            println("NON-EMERGENCY PATIENT NUMBER ==> $number_g")
            println("EMERGENCY PATIENT NUMBER ==> $number_r")
            println("HOSPITAL NUMBER ==> $number_h")
            println("NON-EMERGENCY PATIENT WEIGHTS ==> $WG")
            println("EMERGENCY PATIENT WEIGHTS ==> $WR")
            println("BUDGET ==> $budget")

            # ### Mathematical programming ###
            # formulation_version = "v0929"
            # run_math(problem, 
            #          coordi,
            #          cutting_num, 
            #          number_f, 
            #          number_g,
            #          number_r,
            #          number_h,
            #          WG,
            #          WR,
            #          service,
            #          patient,
            #          demand,
            #          distance_matrix,
            #          Travel,
            #          set_and_parameter[1],
            #          set_and_parameter[2],
            #          set_and_parameter[3],
            #          set_and_parameter[4],
            #          set_and_parameter[5],
            #          set_and_parameter[6],
            #          set_and_parameter[7],
            #          set_and_parameter[8],
            #          set_and_parameter[9],
            #          capa_room, 
            #          capa_man,
            #          set_and_parameter[10],
            #          set_and_parameter[11],
            #          set_and_parameter[12],
            #          set_and_parameter[13],
            #          set_and_parameter[14],
            #          set_and_parameter[15],
            #          budget,
            #          instance_str,
            #          formulation_version)


            ### Heuristic algorithm ###
            heuristic_version = "v0929_clustering_route_split_vehicle_add"
            run_GA_clustering_route_split_vehicle_add(problem, 
                                                    cutting_num,
                                                    number_f, 
                                                    number_g,
                                                    number_r,
                                                    number_h,
                                                    WG,
                                                    WR,
                                                    service,
                                                    patient,
                                                    demand,
                                                    coordi,
                                                    distance_matrix,
                                                    Travel,
                                                    set_and_parameter[1],
                                                    set_and_parameter[2],
                                                    set_and_parameter[3],
                                                    set_and_parameter[4],
                                                    set_and_parameter[5],
                                                    set_and_parameter[7],
                                                    set_and_parameter[8],
                                                    set_and_parameter[9],
                                                    capa_room, 
                                                    capa_man,
                                                    set_and_parameter[10],
                                                    set_and_parameter[11],
                                                    set_and_parameter[12],
                                                    set_and_parameter[13],
                                                    set_and_parameter[14],
                                                    budget,
                                                    instance_str,
                                                    heuristic_version,
                                                    true)

        end
    end
end

prodhon_experiment()