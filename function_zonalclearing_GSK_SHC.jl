m = Model(Gurobi.Optimizer)

function market_clearing!(a::Model)
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()
    m.ext[:objective] = Dict()

    # Extract sets
    N = a.ext[:sets][:N]
    G = a.ext[:sets][:G]
    Z = a.ext[:sets][:Z]
    L = a.ext[:sets][:L]
    L_DC = a.ext[:sets][:L_DC]
    T = a.ext[:sets][:T]

    # Extract parameters
    GENCAP = a.ext[:parameters][:GENCAP]
    MC = a.ext[:parameters][:MC]
    DEM = a.ext[:parameters][:DEM]
    NPTDF = a.ext[:parameters][:NPTDF]
    n_in_z = a.ext[:parameters][:n_in_z]
    g_in_n = a.ext[:parameters][:g_in_n]
    TC = a.ext[:parameters][:TC]
    TC_DC = a.ext[:parameters][:TC_DC]
    inc = a.ext[:parameters][:inc] 
    inc_dc = a.ext[:parameters][:inc_dc]
    RES = a.ext[:parameters][:RES]
    CC = a.ext[:parameters][:CC]

    # Create variables
    v = m.ext[:variables][:v] = @variable(m, [g in G, t in T], lower_bound = 0, upper_bound = 1, base_name = "dispatch")
    curt = m.ext[:variables][:curt] = @variable(m, [n in N, t in T], lower_bound = 0, base_name = "curtailment")
    p = m.ext[:variables][:p] = @variable(m, [n in N, t in T], base_name = "position")
    F = m.ext[:variables][:F] = @variable(m, [l in L, t in T], base_name = "flows")
    F_dc = m.ext[:variables][:F_dc] = @variable(m, [l_dc in L_DC, t in T], base_name = "DC flows")
    
    #Test if parameters/variables are correctly indexed
    # for t in 1:min(2, length(T))  # Adjust accordingly if T is not directly indexable
    #     println("Accessing F_dc[1, $t]: ", F_dc[1,t])
    #     if length(L_DC) >= 2
    #         println("Accessing F_dc[2, $t]: ", F_dc[2,t])
    #     end
    # end
    # for l_dc in L_DC
    #     println("Accessing TC_DC[$l_dc]: ", TC_DC[l_dc])
    # end

    # Create expressions
    GC_per_time_and_zone = m.ext[:expressions][:GC_per_zone] = @expression(m, [z in Z, t in T], sum(sum(MC[g]*GENCAP[g]*v[g,t] for g in g_in_n[n]) + CC*curt[n,t] for n in n_in_z[z]) )
    p_zonal = m.ext[:expressions][:p_zonal] = @expression(m, [z in Z, t in T], sum(sum(GENCAP[g] * v[g,t] for g in g_in_n[n]) + (RES[t,findfirst(N .== n)]-curt[n,t]) - DEM[t,findfirst(N .== n)] for n in n_in_z[z]))            
    
    # Objective
    GC = m.ext[:objective][:GC] = @objective(m, Min, sum(sum(MC[g]*GENCAP[g]*v[g,t] for g in g_in_n[n]) +CC*curt[n,t] for t in T for n in N)) #Minimize the total generation cost (and curtailment costs)


    # Constraints 
    m.ext[:constraints][:con1] = @constraint(m, con1[n in N, t in T], sum(GENCAP[g] * v[g,t] for g in g_in_n[n]) + (RES[t,findfirst(N .== n)]-curt[n,t]) - DEM[t,findfirst(N .== n)] - p[n,t] == 0 ) #Constraint 1d: Power Balance constraint: All Gen cap (for every g, t) + RES (minus curtailed RES) minus the Demand for every hour and node is equal to the net positions for all g and t (or minus p = zero here) 
    m.ext[:constraints][:con3] = @constraint(m, con3[l in L, t in T], F[l,t] == sum(NPTDF[findfirst(N .== n),l]*(sum(GENCAP[g]*(v[g,t]) for g in g_in_n[n]) + (RES[t,findfirst(N .== n)]-curt[n,t])-DEM[t,findfirst(N .== n)] - sum(F_dc[l_dc,t]*inc_dc[l_dc,findfirst(N .== n)] for l_dc in L_DC)) for n in N ))
    m.ext[:constraints][:con4] = @constraint(m, con4[l in L, t in T], -TC[l] <= F[l,t] <= TC[l] )                      #Constraint stating the Flow on the AC line must be between the max Transfer Capacity of that line
    m.ext[:constraints][:con5] = @constraint(m, con5[l_dc in L_DC, t in T], -TC_DC[l_dc] <= F_dc[l_dc,t] <= TC_DC[l_dc] ) #Constraint stating the Flow on the DC line must be between the max Transfer Capacity of that DC line
    m.ext[:constraints][:con6] = @constraint(m, con6[n in N, t in T], curt[n,t] <= RES[t,findfirst(N .== n)] )         #Constraint stating that curtailed power cannot be larger than produced renewable power

    #enter laws of Kirchoff for the curtailment of the OWF in the OBZ. Specified for the node containing the OWF. For a radial connected wind farm use this:
    # read_data = a.ext[:loops][:read_data]
    # if read_data in ["3-node_HM_AC", "3-node_HM_DC","3-node_OBZ_AC","3-node_OBZ_DC"]
    #     m.ext[:constraints][:con9] = @constraint(m, con9[t in T], -F_dc[1,t]-(RES[t,findfirst(N .== 3)]-curt[3,t]) == 0)
    # elseif read_data in ["4-node_OBZ_AC", "4-node_OBZ_DC", "4-node_HM_AC", "4-node_HM_DC"]
    #     m.ext[:constraints][:con9] = @constraint(m, con9[t in T], -F_dc[1,t]-(RES[t,findfirst(N .== 4)]-curt[4,t]) == 0)
    # elseif read_data in ["5-node_HM_AC", "5-node_HM_DC","5-node_OBZ_AC", "5-node_HM_DC"] 
    #     m.ext[:constraints][:con9] = @constraint(m, con9[t in T], -F_dc[1,t]-(RES[t,findfirst(N .== 5)]-curt[5,t]) == 0)
    # end    

    #For a hybrid windfarm use this. Make sure to correctly index the DC lines and the OWF (Only Cross-border DC!)
    read_data = a.ext[:loops][:read_data]
    if read_data in ["3-node_HM_DC","3-node_OBZ_DC"]
        m.ext[:constraints][:con7] = @constraint(m, con7[t in T], -F_dc[1,t]-(RES[t,findfirst(N .== 3)]-curt[3,t]) + F_dc[2,t] == 0)
    elseif read_data in ["4-node_OBZ_DC", "4-node_HM_DC"]
        m.ext[:constraints][:con7] = @constraint(m, con7[t in T], -F_dc[1,t]-(RES[t,findfirst(N .== 4)]-curt[4,t]) + F_dc[2,t] == 0)
    elseif read_data in ["5-node_HM_DC", "5-node_HM_DC"] 
        m.ext[:constraints][:con7] = @constraint(m, con7[t in T], -F_dc[1,t]-(RES[t,findfirst(N .== 5)]-curt[5,t])+ F_dc[2,t] == 0)
    end    

    #OLD Kirchoff law constraints: The DC flow from the shore to the OWF, minus the (RES production -curtailment of OWF) equals the two flows going to the other two wind farms  
    # m.ext[:constraints][:con10] = @constraint(m, con10[t in T], -F_dc[2,t]-(RES[t,findfirst(N .== 119)]-curt[119,t]) +F_dc[8,t] +F_dc[7,t] == 0)
    # m.ext[:constraints][:con11] = @constraint(m, con11[t in T], -F_dc[5,t]-(RES[t,findfirst(N .== 120)]-curt[120,t]) -F_dc[7,t] +F_dc[9,t] == 0)
    # m.ext[:constraints][:con12] = @constraint(m, con12[t in T], -F_dc[6,t]-(RES[t,findfirst(N .== 124)]-curt[124,t]) -F_dc[8,t] -F_dc[9,t] == 0)
    # m.ext[:constraints][:con13] = @constraint(m, con13[t in T], -F_dc[3,t]-(RES[t,findfirst(N .== 121)]-curt[121,t]) +F_dc[10,t] == 0)
    # m.ext[:constraints][:con14] = @constraint(m, con14[t in T], -F_dc[4,t]-(RES[t,findfirst(N .== 123)]-curt[123,t]) -F_dc[10,t] == 0)
    

    return m
end

# Build your model
market_clearing!(a)
status = optimize!(m)

v_bc = a.ext[:parameters][:v_bc] = value.(m.ext[:variables][:v])
curt_bc = a.ext[:parameters][:curt_bc] = value.(m.ext[:variables][:curt])
p_bc = a.ext[:parameters][:p_bc] = value.(m.ext[:expressions][:p_zonal])
F_bc = a.ext[:parameters][:F_bc] = value.(m.ext[:variables][:F])

function get_GSK!(a::Model)
    a.ext[:variables] = Dict()
    a.ext[:expressions] = Dict()
    a.ext[:constraints] = Dict()
    a.ext[:objective] = Dict()

    # Extract sets
    N = a.ext[:sets][:N]
    G = a.ext[:sets][:G]
    Z = a.ext[:sets][:Z]
    L = a.ext[:sets][:L]
    T = a.ext[:sets][:T]

    # Extract parameters
    GENCAP = a.ext[:parameters][:GENCAP]
    n_in_z = a.ext[:parameters][:n_in_z]
    g_in_n = a.ext[:parameters][:g_in_n]
    all_g_in_n = a.ext[:parameters][:all_g_in_n]

    gsk = zeros(length(N),length(Z))
    for n in N
        zone_temp = df_node.Zone[df_node[:,:BusID].==n][1]
        nodes_in_zone = n_in_z[zone_temp]
        if n in nodes_in_zone
            pmax_in_zone = sum(df_gen.Pmax[[x in nodes_in_zone for x in df_gen[:,:OnBus]]])
            pmax_at_node = sum(df_gen.Pmax[[x == n for x in df_gen[:,:OnBus]]])
            gsk[findfirst(N .== n),zone_temp] = pmax_at_node/pmax_in_zone
        end
    end
    
    a.ext[:parameters][:GSK] = gsk

    return a 
end

# Build your model
get_GSK!(a)

ZPTDF, RAM_pos, RAM_neg = zeros(length(Z), length(L)), zeros(length(L), length(T)), zeros(length(L), length(T))
for z in Z
    for l in L
        ZPTDF[z,l] = sum(NPTDF[findfirst(N .== n),l]*a.ext[:parameters][:GSK][findfirst(N .== n),z] for n in N )
    end
end 
for t in T
    for l in L
        RAM_pos[l,t] = max(TC[l] - (F_bc[l,t] - sum(ZPTDF[z,l] * p_bc[z,t] for z in Z )), MinRAM*TC[l])
        RAM_neg[l,t] = max(TC[l] + (F_bc[l,t] - sum(ZPTDF[z,l] * p_bc[z,t] for z in Z )), MinRAM*TC[l])
    end
end
a.ext[:parameters][:ZPTDF] = ZPTDF
a.ext[:parameters][:RAM_pos] = RAM_pos
a.ext[:parameters][:RAM_neg] = RAM_neg

function market_clearing!(a::Model)
    a.ext[:variables] = Dict()
    a.ext[:expressions] = Dict()
    a.ext[:constraints] = Dict()
    a.ext[:objective] = Dict()

    # Extract sets
    N = a.ext[:sets][:N]
    G = a.ext[:sets][:G]
    Z = a.ext[:sets][:Z]
    L = a.ext[:sets][:L]
    L_DC = a.ext[:sets][:L_DC]
    T = a.ext[:sets][:T]

    # Extract parameters
    GENCAP = a.ext[:parameters][:GENCAP]
    MC = a.ext[:parameters][:MC]
    DEM = a.ext[:parameters][:DEM]
    NPTDF = a.ext[:parameters][:NPTDF]
    n_in_z = a.ext[:parameters][:n_in_z]
    g_in_n = a.ext[:parameters][:g_in_n]
    TC = a.ext[:parameters][:TC]
    TC_DC = a.ext[:parameters][:TC_DC]
    v_bc = a.ext[:parameters][:v_bc]
    p_bc = a.ext[:parameters][:p_bc]
    F_bc = a.ext[:parameters][:F_bc]
    GSK = a.ext[:parameters][:GSK]
    inc_dc = a.ext[:parameters][:inc_dc]
    i = a.ext[:parameters][:i]
    i_DC = a.ext[:parameters][:i_DC]
    RES = a.ext[:parameters][:RES]
    CC = a.ext[:parameters][:CC]
    ZPTDF = a.ext[:parameters][:ZPTDF] 
    RAM_pos = a.ext[:parameters][:RAM_pos] 
    RAM_neg = a.ext[:parameters][:RAM_neg] 
    
    # NTC_12, NTC_23, NTC_34, NTC_15, NTC_25, NTC_35, NTC_45 = a.ext[:parameters][:NTC_12], a.ext[:parameters][:NTC_23], a.ext[:parameters][:NTC_34], a.ext[:parameters][:NTC_15], a.ext[:parameters][:NTC_25], a.ext[:parameters][:NTC_35], a.ext[:parameters][:NTC_45] = TC_DC[11], TC_DC[12], TC_DC[13], TC_DC[2], TC_DC[1]+TC_DC[6], TC_DC[4]+TC_DC[5], TC_DC[3]
   
    read_data = a.ext[:loops][:read_data]
    if read_data in ["4-node_HM_DC","4-node_OBZ_DC"]
        NTC_24 = a.ext[:parameters][:NTC_24] = TC_DC[1]
        NTC_34 = a.ext[:parameters][:NTC_34] = TC_DC[2]

    elseif read_data in ["5-node_HM_DC","5-node_OBZ_DC"]
        NTC_35 = a.ext[:parameters][:NTC_35] = TC_DC[1]
        NTC_45 = a.ext[:parameters][:NTC_45] = TC_DC[2]
    end
   
    
    # # Function to dynamically access and potentially use NTC parameters within the model
    # function access_ntc_parameters(a::Model)
    #     if haskey(a.ext[:parameters], :NTC)
    #       for (ntc_key, ntc_value) in a.ext[:parameters][:NTC]
    #             println("$ntc_key: $ntc_value")
    #         end
    #     else
    #         println("No NTC parameters found in the dictionary.")
    #  end
    # end
    # access_ntc_parameters(a)


    # Create variables: Note that variables are slightly differtent to distinghuis the flow-based parameters from the other parameters (e.g. f_dc in model m and now f_DC in model a)
    v = a.ext[:variables][:v] = @variable(a, [g in G, t in T], lower_bound = 0, upper_bound = 1, base_name = "dispatch")
    curt = a.ext[:variables][:curt] = @variable(a, [n in N, t in T], lower_bound = 0, base_name = "curtailment")
    p = a.ext[:variables][:p] = @variable(a, [z in Z, t in T], base_name = "total position")
    p_FB = a.ext[:variables][:p_FB] = @variable(a, [z in Z, t in T], base_name = "flow-based position")
    F_FBMC = a.ext[:variables][:F_FBMC] = @variable(a, [l in L, t in T], base_name = "commercial flow")
    F_DC = a.ext[:variables][:F_DC] = @variable(a, [l_dc in L_DC, t in T], base_name = "commercial DC flow") 

    # Create expressions
    F = a.ext[:expressions][:F] = @expression(a, [l in L, t in T], sum(NPTDF[findfirst(N .== n),l]*(sum(GENCAP[g]*(v[g,t]) for g in g_in_n[n]) + (RES[t,findfirst(N .== n)]-curt[n,t])-DEM[t,findfirst(N .== n)] - sum(F_DC[l_dc,t]*inc_dc[l_dc,findfirst(N .== n)] for l_dc in L_DC)) for n in N ))
    ZPTDF = a.ext[:expressions][:ZPTDF] = @expression(a, [z in Z, l in L], sum(NPTDF[findfirst(N .== n),l]*GSK[findfirst(N .== n),z] for n in N )    )
    RAM_pos = a.ext[:expressions][:RAM_pos] = @expression(a, [l in L, t in T], TC[l] - (F_bc[l,t] - sum(ZPTDF[z,l] * p_bc[z,t] for z in Z )) )
    RAM_neg = a.ext[:expressions][:RAM_neg] = @expression(a, [l in L, t in T], TC[l] + (F_bc[l,t] - sum(ZPTDF[z,l] * p_bc[z,t] for z in Z )) )
    # RAM_pos = a.ext[:expressions][:RAM_pos] = @expression(a, [l in L, t in T], TC[l] )
    # RAM_neg = a.ext[:expressions][:RAM_neg] = @expression(a, [l in L, t in T], TC[l] )
    GC_per_time_and_zone = a.ext[:expressions][:GC_per_zone] = @expression(a, [z in Z, t in T], sum(sum(MC[g]*GENCAP[g]*v[g,t] for g in g_in_n[n]) + CC*curt[n,t] for n in n_in_z[z]) )
    
    # Objective
    GC = a.ext[:objective][:GC] = @objective(a, Min, sum(sum(MC[g]*GENCAP[g]*v[g,t] for g in g_in_n[n]) +CC*curt[n,t] for t in T for n in N))

    # Constraints 
    a.ext[:constraints][:con1] = @constraint(a, con1[z in Z, t in T], sum(sum(GENCAP[g] * v[g,t] for g in g_in_n[n]) + (RES[t,findfirst(N .== n)]-curt[n,t]) - DEM[t,findfirst(N .== n)] for n in n_in_z[z]) - p[z,t] == 0 )
    a.ext[:constraints][:con2] = @constraint(a, con2[z in Z, t in T], sum(F_FBMC[l,t]*i[l,z] for l in L) + sum(F_DC[l_dc,t]*i_DC[l_dc,z] for l_dc in L_DC) - p[z,t] == 0 )
    a.ext[:constraints][:con3] = @constraint(a, con3[l in L, t in T], F_FBMC[l,t] == sum(ZPTDF[z,l]*p_FB[z,t] for z in Z))
    a.ext[:constraints][:con4] = @constraint(a, con4[z in Z, t in T], p_FB[z,t] == p[z,t] - sum(F_DC[l_dc,t]*i_DC[l_dc,z] for l_dc in L_DC))
    a.ext[:constraints][:con5] = @constraint(a, con5[n in N, t in T], curt[n,t] <= RES[t,findfirst(N .== n)] )
    a.ext[:constraints][:con6] = @constraint(a, con6[l in L, t in T], -RAM_neg[l,t] <= F_FBMC[l,t] <= RAM_pos[l,t] )

    #enter laws of Kirchoff: Same as above, only attributed to model a instead of model m
    read_data = a.ext[:loops][:read_data]
    if read_data in ["3-node_HM_DC","3-node_OBZ_DC"]
        m.ext[:constraints][:con7] = @constraint(a, con7[t in T], -F_DC[1,t]-(RES[t,findfirst(N .== 3)]-curt[3,t]) - F_DC[2,t] == 0)
    elseif read_data in ["4-node_OBZ_DC", "4-node_HM_DC"]
        m.ext[:constraints][:con7] = @constraint(a, con7[t in T], -F_DC[1,t]-(RES[t,findfirst(N .== 4)]-curt[4,t]) - F_DC[2,t] == 0)
    elseif read_data in ["5-node_HM_DC", "5-node_HM_DC"] 
        m.ext[:constraints][:con7] = @constraint(a, con7[t in T], -F_DC[1,t]-(RES[t,findfirst(N .== 5)]-curt[5,t]) - F_DC[2,t] == 0)
    end   

    
    # #enter laws of Kirchoff
    # a.ext[:constraints][:con9] = @constraint(a, con9[t in T], -F_DC[1,t]-(RES[t,findfirst(N .== 122)]-curt[122,t]) == 0)
    # a.ext[:constraints][:con10] = @constraint(a, con10[t in T], -F_DC[2,t]-(RES[t,findfirst(N .== 119)]-curt[119,t]) +F_DC[8,t] +F_DC[7,t] == 0)
    # a.ext[:constraints][:con11] = @constraint(a, con11[t in T], -F_DC[5,t]-(RES[t,findfirst(N .== 120)]-curt[120,t]) -F_DC[7,t] +F_DC[9,t] == 0)
    # a.ext[:constraints][:con12] = @constraint(a, con12[t in T], -F_DC[6,t]-(RES[t,findfirst(N .== 124)]-curt[124,t]) -F_DC[8,t] -F_DC[9,t] == 0)
    # a.ext[:constraints][:con13] = @constraint(a, con13[t in T], -F_DC[3,t]-(RES[t,findfirst(N .== 121)]-curt[121,t]) +F_DC[10,t] == 0)
    # a.ext[:constraints][:con14] = @constraint(a, con14[t in T], -F_DC[4,t]-(RES[t,findfirst(N .== 123)]-curt[123,t]) -F_DC[10,t] == 0)

    if read_data in ["4-node_HM_DC","4-node_OBZ_DC"]
        a.ext[:constraints][:con8] = @constraint(a, con8[t in T], -NTC_24 <= F_DC[1,t] <= NTC_24) #limit flow on cross-border line
        a.ext[:constraints][:con9] = @constraint(a, con9[t in T], -NTC_34 <= F_DC[2,t] <= NTC_34) #limit flow on cross-border line
    
    elseif read_data in ["5-node_HM_DC","5-node_OBZ_DC"]
        a.ext[:constraints][:con8] = @constraint(a, con8[t in T], -NTC_35 <= F_DC[1,t] <= NTC_35) #limit flow on cross-border line
        a.ext[:constraints][:con9] = @constraint(a, con9[t in T], -NTC_45 <= F_DC[2,t] <= NTC_45) #limit flow on cross-border line
    end

   #NTC values OLD
#    a.ext[:constraints][:con15] = @constraint(a, con15[t in T], -NTC_12 <= F_DC[11,t] <= NTC_12) #limit flow on cross-border line
#    a.ext[:constraints][:con16] = @constraint(a, con16[t in T], -NTC_23 <= -F_DC[12,t] <= NTC_23) #limit flow on cross-border line
#    a.ext[:constraints][:con17] = @constraint(a, con17[t in T], -NTC_34 <= F_DC[13,t] <= NTC_34) #limit flow on cross-border line
#    a.ext[:constraints][:con18] = @constraint(a, con18[t in T], -NTC_15 <= F_DC[2,t] <= NTC_15) #limit flow on cross-border line
#    a.ext[:constraints][:con19] = @constraint(a, con19[t in T], -NTC_25 <= F_DC[1,t]+F_DC[6,t] <= NTC_25) #limit flow on cross-border line
#    a.ext[:constraints][:con20] = @constraint(a, con20[t in T], -NTC_35 <= F_DC[4,t]+F_DC[5,t] <= NTC_35) #limit flow on cross-border line
#    a.ext[:constraints][:con21] = @constraint(a, con21[t in T], -NTC_45 <= F_DC[3,t] <= NTC_45) #limit flow on cross-border line
 
   # NTC values
    # a.ext[:constraints][:con15] = @constraint(a, con15[t in T], -NTC_12 <= F_DC[11,t]+F_DC[8,t] <= NTC_12) #limit flow on cross-border line
    # a.ext[:constraints][:con16] = @constraint(a, con16[t in T], -NTC_23 <= -F_DC[12,t]-F_DC[9,t] <= NTC_23) #limit flow on cross-border line
    # a.ext[:constraints][:con17] = @constraint(a, con17[t in T], -NTC_34 <= F_DC[13,t]-F_DC[10,t] <= NTC_34) #limit flow on cross-border line
    # a.ext[:constraints][:con18] = @constraint(a, con18[t in T], -NTC_13 <= F_DC[7,t] <= NTC_13) #limit flow on cross-border line



    return a
end

# Build your model
market_clearing!(a)
@info "Market clearing model built."
status = optimize!(a)
@info "Market clearing model optimised."

GC_DA = a.ext[:parameters][:GC_DA] = value.(a.ext[:expressions][:GC_per_zone])
v_DA = a.ext[:parameters][:v_DA] = value.(a.ext[:variables][:v])
curt_DA = a.ext[:parameters][:curt_DA] = value.(a.ext[:variables][:curt])
p_DA = a.ext[:parameters][:p_DA] = value.(a.ext[:variables][:p])
F_DA = a.ext[:parameters][:F_DA] = value.(a.ext[:expressions][:F])
F_DC_DA = a.ext[:parameters][:F_DC_DA] = value.(a.ext[:variables][:F_DC])
RES_prod = a.ext[:parameters][:RES_prod] = value.(a.ext[:parameters][:RES])



for l in L
    if (value.(a.ext[:variables][:F_FBMC])[l, 1] == a.ext[:parameters][:RAM_pos]) || (value.(a.ext[:variables][:F_FBMC])[l, 1] == a.ext[:parameters][:RAM_neg])
        println("line ", l, " is congested")
    else
        println("AC line ", l, " is not congested")    
    end

end

read_data = a.ext[:loops][:read_data]
if read_data in ["4-node_HM_DC","4-node_OBZ_DC"]
    for l_dc in L_DC
        if (value.(a.ext[:variables][:F_DC])[l_dc, 1] == a.ext[:parameters][:NTC_24]) || (value.(a.ext[:variables][:F_DC])[l_dc, 1] == -a.ext[:parameters][:NTC_34])
            println("DC line ", l_dc, " congested")
        else
            println("DC line ", l_dc, " is not congested.")
        end
    end
elseif read_data in ["5-node_HM_DC","5-node_OBZ_DC"]
    for l_dc in L_DC
        if (value.(a.ext[:variables][:F_DC])[l_dc, 1] == a.ext[:parameters][:NTC_35]) || (value.(a.ext[:variables][:F_DC])[l_dc, 1] == -a.ext[:parameters][:NTC_45])
            println("DC line ", l_dc, " congested")
        else
            println("DC line ", l_dc, " is not congested.")
        end
    end
end
