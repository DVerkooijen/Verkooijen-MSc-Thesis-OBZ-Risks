m = Model(Gurobi.Optimizer)


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
    inc_dc = a.ext[:parameters][:inc_dc]
    i = a.ext[:parameters][:i]
    i_DC = a.ext[:parameters][:i_DC]
    RES = a.ext[:parameters][:RES]
    CC = a.ext[:parameters][:CC]

    # Create variables
    v = a.ext[:variables][:v] = @variable(a, [g in G, t in T], lower_bound = 0, upper_bound = 1, base_name = "dispatch")            #Constraint 1b (Kenis, 2023)
    vbar = a.ext[:variables][:vbar] = @variable(a, [g in G, t in T], lower_bound = 0, upper_bound = 1, base_name = "dispatch bar")  #constraint 2c (Kenis, 2023)
    curt = a.ext[:variables][:curt] = @variable(a, [n in N, t in T], lower_bound = 0, base_name = "curtailment")                    #constraint 1c  lower bound (Kenis, 2023)
    curtbar = a.ext[:variables][:curtbar] = @variable(a, [n in N, t in T], lower_bound = 0, base_name = "curtailmentbar")           #constraint 2b (Kenis, 2023)
    p = a.ext[:variables][:p] = @variable(a, [z in Z, t in T], base_name = "total position")                                        #p_z (Kenis, 2023)
    F_DC = a.ext[:variables][:F_DC] = @variable(a, [l_dc in L_DC, t in T], base_name = "commercial DC flow")                        #small f,l^DC (Kenis, 2023)
    F_DCbar = a.ext[:variables][:F_DCbar] = @variable(a, [l_dc in L_DC, t in T], base_name = "commercial DC flow bar")              #small f,l^DC' (Kenis, 2023)


    # Create expressions
    F = a.ext[:expressions][:F] = @expression(a, [l in L, t in T], sum(NPTDF[findfirst(N .== n),l]*(sum(GENCAP[g]*(v[g,t]) for g in g_in_n[n]) + (RES[t,findfirst(N .== n)]-curt[n,t])-DEM[t,findfirst(N .== n)] - sum(F_DC[l_dc,t]*inc_dc[l_dc,findfirst(N .== n)] for l_dc in L_DC)) for n in N )) 
    Fbar = a.ext[:expressions][:Fbar] = @expression(a, [l in L, t in T], sum(NPTDF[findfirst(N .== n),l]*(sum(GENCAP[g]*(vbar[g,t]) for g in g_in_n[n]) + (RES[t,findfirst(N .== n)]-curtbar[n,t])-DEM[t,findfirst(N .== n)] - sum(F_DCbar[l_dc,t]*inc_dc[l_dc,findfirst(N .== n)] for l_dc in L_DC)) for n in N )) #Part of constraint 2e (Kenis, 2023)
    
    GC_per_time_and_zone = a.ext[:expressions][:GC_per_zone] = @expression(a, [z in Z, t in T], sum(sum(MC[g]*GENCAP[g]*v[g,t] for g in g_in_n[n]) + CC*curt[n,t] for n in n_in_z[z]) )

    # Objective
    GC = a.ext[:objective][:GC] = @objective(a, Min, sum(sum(MC[g]*GENCAP[g]*v[g,t] for g in g_in_n[n]) +CC*curt[n,t] for t in T for n in N)) #Objective function 1a (Kenis, 2023)

    # Constraints 
    a.ext[:constraints][:con1] = @constraint(a, con1[z in Z, t in T], sum(sum(GENCAP[g] * v[g,t] for g in g_in_n[n]) + (RES[t,findfirst(N .== n)]-curt[n,t]) - DEM[t,findfirst(N .== n)] for n in n_in_z[z]) - p[z,t] == 0 )    #Constraint 1d: Power balance (Kenis, 2023)
   
    a.ext[:constraints][:con2] = @constraint(a, con2[z in Z, t in T], sum(F[l,t]*i[l,z] for l in L) + sum(F_DC[l_dc,t]*i_DC[l_dc,z] for l_dc in L_DC) - p[z,t] == 0 )  #Constraint 1e (Kenis, 2023)
    a.ext[:constraints][:con3] = @constraint(a, con3[n in N, t in T], curt[n,t] <= RES[t,findfirst(N .== n)] )  #constraint 1c upper bound (Kenis, 2023)
  
    a.ext[:constraints][:con4] = @constraint(a, con4[z in Z, t in T], sum(sum(GENCAP[g] * vbar[g,t] for g in g_in_n[n]) + (RES[t,findfirst(N .== n)]-curtbar[n,t]) - DEM[t,findfirst(N .== n)] for n in n_in_z[z]) - p[z,t] == 0 )      #Constraint 2c: Power balance projected (Kenis, 2023)
    a.ext[:constraints][:con5] = @constraint(a, con5[z in Z, t in T], sum(Fbar[l,t]*i[l,z] for l in L) + sum(F_DCbar[l_dc,t]*i_DC[l_dc,z] for l_dc in L_DC) - p[z,t] == 0 ) #Constraint 2d (Kenis, 2023)
    a.ext[:constraints][:con6] = @constraint(a, con6[l in L, t in T], -TC[l] <= Fbar[l,t] <= TC[l]) #Coinstaint 2e (Kenis, 2023)
    a.ext[:constraints][:con7] = @constraint(a, con7[l_dc in L_DC, t in T], -TC_DC[l_dc] <= F_DCbar[l_dc,t] <= TC_DC[l_dc]) #constaint 2f (Kenis, 2023)
    a.ext[:constraints][:con8] = @constraint(a, con8[n in N, t in T], curtbar[n,t] <= RES[t,findfirst(N .== n)] )   #Constraint 2b upper bound (Kenis, 2023)

    #enter laws of Kirchoff
    a.ext[:constraints][:con9] = @constraint(a, con9[t in T], -F_DCbar[1,t]-(RES[t,findfirst(N .== 5)]-curtbar[5,t]) - F_DCbar[2,t] == 0)       #Constraint 2g (Kenis, 2023): nodal power balance for nodes only connected to DC lines (The OWFs in this case!)
    
    #enter laws of Kirchoff
    a.ext[:constraints][:con10] = @constraint(a, con10[t in T], -F_DC[1,t]-(RES[t,findfirst(N .== 5)]-curt[5,t]) - F_DC[2,t] == 0)  #


    # a.ext[:constraints][:con15] = @constraint(a, con15[t in T], -F_DC[1,t]-(RES[t,findfirst(N .== 122)]-curt[122,t]) == 0)
    # a.ext[:constraints][:con16] = @constraint(a, con16[t in T], -F_DC[2,t]-(RES[t,findfirst(N .== 119)]-curt[119,t]) +F_DC[8,t] +F_DC[7,t] == 0)
    # a.ext[:constraints][:con17] = @constraint(a, con17[t in T], -F_DC[5,t]-(RES[t,findfirst(N .== 120)]-curt[120,t]) -F_DC[7,t] +F_DC[9,t] == 0)
    # a.ext[:constraints][:con18] = @constraint(a, con18[t in T], -F_DC[6,t]-(RES[t,findfirst(N .== 124)]-curt[124,t]) -F_DC[8,t] -F_DC[9,t] == 0)
    # a.ext[:constraints][:con19] = @constraint(a, con19[t in T], -F_DC[3,t]-(RES[t,findfirst(N .== 121)]-curt[121,t]) +F_DC[10,t] == 0)
    # a.ext[:constraints][:con20] = @constraint(a, con20[t in T], -F_DC[4,t]-(RES[t,findfirst(N .== 123)]-curt[123,t]) -F_DC[10,t] == 0)

    # a.ext[:constraints][:con9] = @constraint(a, con9[t in T], -F_DCbar[1,t]-(RES[t,findfirst(N .== 122)]-curtbar[122,t]) == 0)
    # a.ext[:constraints][:con10] = @constraint(a, con10[t in T], -F_DCbar[2,t]-(RES[t,findfirst(N .== 119)]-curtbar[119,t]) +F_DCbar[8,t] +F_DCbar[7,t] == 0)
    # a.ext[:constraints][:con11] = @constraint(a, con11[t in T], -F_DCbar[5,t]-(RES[t,findfirst(N .== 120)]-curtbar[120,t]) -F_DCbar[7,t] +F_DCbar[9,t] == 0)
    # a.ext[:constraints][:con12] = @constraint(a, con12[t in T], -F_DCbar[6,t]-(RES[t,findfirst(N .== 124)]-curtbar[124,t]) -F_DCbar[8,t] -F_DCbar[9,t] == 0)
    # a.ext[:constraints][:con13] = @constraint(a, con13[t in T], -F_DCbar[3,t]-(RES[t,findfirst(N .== 121)]-curtbar[121,t]) +F_DCbar[10,t] == 0)
    # a.ext[:constraints][:con14] = @constraint(a, con14[t in T], -F_DCbar[4,t]-(RES[t,findfirst(N .== 123)]-curtbar[123,t]) -F_DCbar[10,t] == 0)

    return a
end

# Build your model
market_clearing!(a)
write_to_file(a, "model_check/DA_exact_pre_optimization.lp")

status = optimize!(a)
write_to_file(a, "model_check/DA_exact_post_optimization.lp")

GC_DA = a.ext[:parameters][:GC_DA] = value.(a.ext[:expressions][:GC_per_zone])
v_DA = a.ext[:parameters][:v_DA] = value.(a.ext[:variables][:v])                #variable dispatch (v_g)
curt_DA = a.ext[:parameters][:curt_DA] = value.(a.ext[:variables][:curt])       #variable curtialment (c_n)
p_DA = a.ext[:parameters][:p_DA] = value.(a.ext[:variables][:p])                #Net position zone (p_z)
F_DA = a.ext[:parameters][:F_DA] = value.(a.ext[:expressions][:F])              #Flow on AC lines
F_DA_bar = a.ext[:parameters][:F_DA_bar] = value.(a.ext[:expressions][:Fbar])   #Flow on AC lines ' (bar)  
F_DC_DA = a.ext[:parameters][:F_DC_DA] = value.(a.ext[:variables][:F_DC])           #Commercial DC flow (f_l^DC)
F_DC_DA_bar = a.ext[:parameters][:F_DC_DA_bar] = value.(a.ext[:variables][:F_DCbar]) #Commercial DC flow bar (f_l^DC')
RES_prod = a.ext[:parameters][:RES_prod] = value.(a.ext[:parameters][:RES])         #RES production

