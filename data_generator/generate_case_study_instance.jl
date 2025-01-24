include("data_generator.jl")


function save_instance(F, G, R, H, TMC_coord, non_emer_coord, emer_coord, hospital_coord, Service, Patient, Demand, budget, capa_room, capa_man)

    log_dir = "data_generator/case_study_instance/"
    mkpath(log_dir)

    name = "FGRH_$(F)_$(G)_$(R)_$(H).txt"
    file = open(string(log_dir, name), "w")
    
    println(file, F)
    println(file, G)
    println(file, R)
    println(file, H)
    println(file, "")

    shuffle!(non_emer_coord)
    shuffle!(emer_coord)

    for (x, y) in TMC_coord 
        println(file, "$(x)\t$(y)")
    end
    println(file, "")

    for (idx, (x, y)) in enumerate(non_emer_coord)
        if idx > G break end
        println(file, "$(x)\t$(y)")
    end
    println(file, "")

    for (idx, (x, y)) in enumerate(emer_coord)
        if idx > R break end
        println(file, "$(x)\t$(y)")
    end
    println(file, "")

    for (idx, (x, y)) in enumerate(hospital_coord)
        if idx > H break end
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


function generate_instance()

    instance_set = [[25, 100, 200, 20, 500000]]
    # instance_set = [[25, 125, 375, 30, 750000], [25, 250, 250, 20, 625000], [25, 375, 125, 10, 500000],
    #                 [25, 250, 750, 70, 1250000], [25, 500, 500, 50, 1000000], [25, 750, 250, 20, 750000]]
    
    # TMC coordinates
    s = open("data_generator/case_study_instance/coord/TMC_coord.txt") do file
        readlines(file)
    end

    TMC_coord = []
    for coor in s
        lat, long = split(coor, "\t")
        push!(TMC_coord, [lat, long])
    end

    # non-emergency patients coordinates
    s = open("data_generator/case_study_instance/coord/non_emer_coord.txt") do file
        readlines(file)
    end

    non_emer_coord = []
    for coor in s
        lat, long = split(coor, "\t")
        push!(non_emer_coord, [lat, long])
    end
    
    # emergency patients coordinates
    s = open("data_generator/case_study_instance/coord/emer_coord.txt") do file
        readlines(file)
    end

    emer_coord = []
    for coor in s
        lat, long = split(coor, "\t")
        push!(emer_coord, [lat, long])
    end

    # hospital coordinates
    s = open("data_generator/case_study_instance/coord/hospital_coord.txt") do file 
        readlines(file)
    end

    hospital_coord = []
    for coor in s
        lat, long = split(coor, "\t")
        push!(hospital_coord, [lat, long])
    end

    for (TMC_num, green_num, red_num, hospital_num, budget) in instance_set

        cutting_num = TMC_num + green_num + red_num + hospital_num
        N = [i for i in 1:cutting_num+TMC_num]
        F = generate_node(N, TMC_num)
        G = generate_node(difference(N, F), green_num)
        R = generate_node(difference(N, vcat(F, G)), red_num)
        H = generate_node(difference(N, vcat(F, G, R)), hospital_num)

        Service = generate_service(N, F, G, R, H)
        Patient = generate_patient(fill(0, cutting_num), G, R)
        Demand = generate_demand(fill(0, cutting_num), vcat(G, R))
        
        capa_room = generate_h_room(Patient, R, H)
        capa_man = generate_h_man(Demand, R, H)
        
        save_instance(TMC_num, green_num, red_num, hospital_num, TMC_coord, non_emer_coord, emer_coord, hospital_coord, Service, Patient, Demand, budget, capa_room, capa_man)
    end
end


generate_instance()