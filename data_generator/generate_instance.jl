include("data_generator.jl")

function Prodhon_LRP_dataset(folder_path)
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

        number_cus = problem[1]
        number_depo = problem[2]
    
        coordi = Vector{Vector}()
        for (x, y) in problem[3:3+number_cus+number_depo-1]
            push!(coordi, [x,y])
        end
    
        vehicle_capa = problem[3+number_cus+number_depo]
        depot_capa = Vector{Int}()
        for c in problem[4+number_cus+number_depo:4+number_cus*2+number_depo-1]
            push!(depot_capa, c)
        end
    
        demands = Vector{Int}()
        for d in problem[4+number_cus*2+number_depo:4+number_cus*2+number_depo*2-1]
            push!(demands, d)
        end
    
        open_cost_depo = Vector{Int}()
        for c in problem[4+number_cus*2+number_depo*2:4+number_cus*2+number_depo*3-1]
            push!(open_cost_depo, c)
        end
    
        open_cost_route = problem[4+number_cus*2+number_depo*3]
    
        push!(instance_set, [file, number_cus, number_depo, coordi])
    end

    return instance_set
end

function save_instance(instance_str, TMC_num, green_num, red_num, hospital_num, coordinates, Service, Patient, Demand, budget, capa_room, capa_man)

    log_dir = "data_generator/TMCLRP_Instance/"
    mkpath(log_dir)

    file = open(string(log_dir, instance_str), "w")
    
    println(file, TMC_num)
    println(file, green_num)
    println(file, red_num)
    println(file, hospital_num)
    println(file, "")

    shuffle!(coordinates)

    for (x, y) in coordinates[1:TMC_num] 
        println(file, "$(x)\t$(y)")
    end
    println(file, "")

    for (x, y) in coordinates[1+TMC_num:TMC_num+green_num] 
        println(file, "$(x)\t$(y)")
    end
    println(file, "")

    for (x, y) in coordinates[1+TMC_num+green_num:TMC_num+green_num+red_num] 
        println(file, "$(x)\t$(y)")
    end
    println(file, "")

    for (x, y) in coordinates[1+TMC_num+green_num+red_num:TMC_num+green_num+red_num+hospital_num] 
        println(file, "$(x)\t$(y)")
    end
    println(file, "")

    for value in Service
        println(file, value)
    end
    println(file, "")

    for value in Patient
        println(file, value)
    end
    println(file, "")

    for value in Demand
        println(file, value)
    end
    println(file, "")

    println(file, budget)
    println(file, "")

    for value in capa_room
        println(file, value)
    end
    println(file, "")

    for value in capa_man
        println(file, value)
    end
    println(file, "")

    close(file)

end


function generate_instance(red_patient_rate, hospital_rate)
    folder_path = "data_generator/Prodhon_LRP" 
    instance_set = Prodhon_LRP_dataset(folder_path)
    
    for (problem, patient_num, TMC_num, coordi) in instance_set

        green_num = round(Int, patient_num*(1-red_patient_rate))
        red_num = patient_num-green_num
        hospital_num = ceil(Int, TMC_num*hospital_rate)
        cutting_num = TMC_num + patient_num + hospital_num

        hospital_coor = generate_hospital(coordi, hospital_num)
        coordinates = vcat(coordi, hospital_coor)

        N = [i for i in 1:cutting_num+TMC_num]
        F = generate_node(N, TMC_num)
        G = generate_node(difference(N, F), green_num)
        R = generate_node(difference(N, vcat(F, G)), red_num)
        H = generate_node(difference(N, vcat(F, G, R)), hospital_num)

        Service = generate_service(N, F, G, R, H)
        Patient = generate_patient(fill(0, cutting_num), G, R)
        Demand = generate_demand(fill(0, cutting_num), vcat(G,R))
        
        budget = 0
        if TMC_num == 5
            if patient_num == 20
                budget = 100000
            elseif patient_num == 50
                budget = 150000
            else
                budget = 200000
            end
        else
            if patient_num == 100
                budget = 250000
            else
                budget = 300000
            end
        end

        capa_room = generate_h_room(Patient, R, H)
        capa_man = generate_h_man(Demand, R, H)

        instance_str = "FGRH_$(TMC_num)_$(green_num)_$(red_num)_$(hospital_num).txt"

        save_instance(instance_str, TMC_num, green_num, red_num, hospital_num, coordinates, Service, Patient, Demand, budget, capa_room, capa_man)
    end
end



red_patient_rate = 0.25
hospital_rate = 0.2

generate_instance(red_patient_rate, hospital_rate)