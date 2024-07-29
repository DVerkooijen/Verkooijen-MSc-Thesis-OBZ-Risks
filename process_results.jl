select_model = a.ext[:loops][:select_model]  
read_data = a.ext[:loops][:read_data]
redispatch = a.ext[:loops][:redispatch]
hydrogen = a.ext[:loops][:hydrogen]
base_folder_results = a.ext[:loops][:base_folder_results]
cap_calc = a.ext[:loops][:cap_calc]

# Optimal values for day-ahead market outcome
GC_DA = a.ext[:parameters][:GC_DA]
GC_opt = a.ext[:parameters][:GC_opt] = zeros(length(a.ext[:sets][:T]))

if redispatch == "yes"
    RDC_perzone_opt = value.(b.ext[:expressions][:RDC_per_zone])
    RDC_opt = a.ext[:parameters][:RDC_opt] = zeros(length(a.ext[:sets][:T]))
    λ_rd_1 = JuMP.dual.(b.ext[:constraints][:con1])
    λ_rd_2 = JuMP.dual.(b.ext[:constraints][:con2])
    λ_rd_3 = JuMP.dual.(b.ext[:constraints][:con3])
    λ_rd_4 = JuMP.dual.(b.ext[:constraints][:con4])
    λ_rd_5 = JuMP.dual.(b.ext[:constraints][:con5])
    λ_rd_6 = JuMP.dual.(b.ext[:constraints][:con6])
    λ_rd_7 = JuMP.dual.(b.ext[:constraints][:con7])
    λ_rd_8 = JuMP.dual.(b.ext[:constraints][:con8])
    if read_data == "5-node_HM" || read_data == "5-node_OBZ"
        λ_rd_9 = JuMP.dual.(b.ext[:constraints][:con9])
    elseif read_data == "Schonheit_HM" || read_data == "Schonheit_OBZ" || read_data == "Schonheit_OBZ_adjusted" || read_data == "Reference_Case" || read_data == "Simple_Hybrid"
        λ_rd_9 = JuMP.dual.(b.ext[:constraints][:con9])
        λ_rd_10 = JuMP.dual.(b.ext[:constraints][:con10])
        λ_rd_11 = JuMP.dual.(b.ext[:constraints][:con11])
        λ_rd_12 = JuMP.dual.(b.ext[:constraints][:con12])
        λ_rd_13 = JuMP.dual.(b.ext[:constraints][:con13])
        if read_data == "Schonheit_HM" || read_data == "Schonheit_OBZ" || read_data == "Schonheit_OBZ_adjusted" || read_data == "Reference_Case"
            λ_rd_14 = JuMP.dual.(b.ext[:constraints][:con14])
        end
    end
    if hydrogen == "yes"
        λ_rd_h2_1 = JuMP.dual.(b.ext[:constraints][:conh2_1])
        λ_rd_h2_2 = JuMP.dual.(b.ext[:constraints][:conh2_2])
    end
end

for t in T
    GC_opt[t] = a.ext[:parameters][:GC_opt][t] = sum(GC_DA[z,t] for z in Z)
    if redispatch == "yes"
        RDC_opt[t] = a.ext[:parameters][:RDC_opt][t] = sum(RDC_perzone_opt[z,t] for z in Z)
    end
end

v_opt = JuMP.value.(a.ext[:variables][:v])
p_opt = JuMP.value.(a.ext[:variables][:p])
curt_opt = a.ext[:parameters][:curt_DA] 
F_AC_opt = a.ext[:parameters][:F_DA] 
F_DC_opt = a.ext[:parameters][:F_DC_DA]
NPTDF_opt = a.ext[:parameters][:NPTDF]



if select_model == "Exact_HC"
    vbar_opt = JuMP.value.(a.ext[:variables][:vbar])
    curtbar_opt = value.(a.ext[:variables][:curtbar])
    F_bar_opt = JuMP.value.(a.ext[:expressions][:Fbar])
    F_DC_opt_bar = a.ext[:parameters][:F_DC_DA_bar]
elseif select_model == "GSK_AHC"
    p_FB_DA = a.ext[:parameters][:p_FB] = value.(a.ext[:variables][:p_FB])
    F_FBMC_DA = a.ext[:parameters][:F_FBMC] = value.(a.ext[:variables][:F_FBMC])

    GC_bc = a.ext[:parameters][:GC_DA] = value.(m.ext[:expressions][:GC_per_zone])
    v_bc = a.ext[:parameters][:v_bc] = value.(m.ext[:variables][:v])
    curt_bc = a.ext[:parameters][:curt_bc] = value.(m.ext[:variables][:curt])
    p_bc = a.ext[:parameters][:p_bc] = value.(m.ext[:expressions][:p_zonal])
    F_bc = a.ext[:parameters][:F_bc] = value.(m.ext[:variables][:F])
    F_dc_bc = a.ext[:parameters][:F_DC_DA] = value.(m.ext[:variables][:F_dc])
    # ZPTDF_opt = a.ext[:parameters][:ZPTDF] = value.(a.ext[:expressions][:ZPTDF])
    # RAM_pos_opt = a.ext[:parameters][:RAM_pos] = value.(a.ext[:expressions][:RAM_pos])
    # RAM_neg_opt = a.ext[:parameters][:RAM_neg] = value.(a.ext[:expressions][:RAM_neg])
    ZPTDF_opt_bc = a.ext[:parameters][:ZPTDF_bc] = value.(a.ext[:parameters][:ZPTDF_bc])
    RAM_pos_opt_bc = a.ext[:parameters][:RAM_pos_bc] = value.(a.ext[:parameters][:RAM_pos_bc])
    RAM_neg_opt_bc = a.ext[:parameters][:RAM_neg_bc] = value.(a.ext[:parameters][:RAM_neg_bc])
    Flow_FBMC_opt = a.ext[:parameters][:Flow_FBMC] = value.(a.ext[:expressions][:Flow_FBMC])
end

if hydrogen == "yes"
    e_DA = a.ext[:parameters][:e_DA]
    e_bc = e_bc = a.ext[:parameters][:e_bc]
end

λ_opt = JuMP.dual.(a.ext[:constraints][:con1])  
λ_opt2 = JuMP.dual.(a.ext[:constraints][:con2])
λ_opt3 = JuMP.dual.(a.ext[:constraints][:con3])
λ_opt4 = JuMP.dual.(a.ext[:constraints][:con4])
λ_opt5 = JuMP.dual.(a.ext[:constraints][:con5])


if read_data == "Schonheit_HM" || read_data == "Schonheit_OBZ" || read_data == "Schonheit_OBZ_adjusted" || read_data == "Reference_Case" || read_data == "Simple_Hybrid"
    λ_opt6 = JuMP.dual.(a.ext[:constraints][:con6])
    λ_opt7 = JuMP.dual.(a.ext[:constraints][:con7])
    λ_opt8 = JuMP.dual.(a.ext[:constraints][:con8])
    λ_opt9 = JuMP.dual.(a.ext[:constraints][:con9])
    λ_opt10 = JuMP.dual.(a.ext[:constraints][:con10])
    λ_opt11 = JuMP.dual.(a.ext[:constraints][:con11])
    λ_opt12 = JuMP.dual.(a.ext[:constraints][:con12])
    λ_opt13 = JuMP.dual.(a.ext[:constraints][:con13])
    λ_opt14 = JuMP.dual.(a.ext[:constraints][:con14])
    λ_opt15 = JuMP.dual.(a.ext[:constraints][:con15])
    λ_opt16 = JuMP.dual.(a.ext[:constraints][:con16])
    if read_data == "Schonheit_HM" || read_data == "Schonheit_OBZ" || read_data == "Schonheit_OBZ_adjusted" || read_data == "Reference_Case"
        λ_opt17 = JuMP.dual.(a.ext[:constraints][:con17])
        λ_opt18 = JuMP.dual.(a.ext[:constraints][:con18])
        λ_opt19 = JuMP.dual.(a.ext[:constraints][:con19])
        λ_opt20 = JuMP.dual.(a.ext[:constraints][:con20])
    end

elseif read_data == "5-node_HM" || read_data == "5-node_OBZ"
    λ_opt6 = JuMP.dual.(a.ext[:constraints][:con6])
    λ_opt7 = JuMP.dual.(a.ext[:constraints][:con7])
    λ_opt8 = JuMP.dual.(a.ext[:constraints][:con8])
end    

if select_model == "GSK_AHC"
    λ_bc_1 = JuMP.dual.(m.ext[:constraints][:con1])
    λ_bc_2 = JuMP.dual.(m.ext[:constraints][:con2])
    λ_bc_3 = JuMP.dual.(m.ext[:constraints][:con3])
    λ_bc_4 = JuMP.dual.(m.ext[:constraints][:con4])
    λ_bc_5 = JuMP.dual.(m.ext[:constraints][:con5])
    λ_bc_6 = JuMP.dual.(m.ext[:constraints][:con6])
    
    if read_data == "Schonheit_HM" || read_data == "Schonheit_OBZ" || read_data == "Schonheit_OBZ_adjusted" || read_data == "Reference_Case" || read_data == "Simple_Hybrid"
        λ_bc_7 = JuMP.dual.(m.ext[:constraints][:con7])
        λ_bc_8 = JuMP.dual.(m.ext[:constraints][:con8])
        λ_bc_9 = JuMP.dual.(m.ext[:constraints][:con9])
        λ_bc_10 = JuMP.dual.(m.ext[:constraints][:con10])
        λ_bc_11 = JuMP.dual.(m.ext[:constraints][:con11])
        λ_bc_12 = JuMP.dual.(m.ext[:constraints][:con12])
        if read_data == "Schonheit_HM" || read_data == "Schonheit_OBZ" || read_data == "Schonheit_OBZ_adjusted" || read_data == "Reference_Case"
            λ_bc_13 = JuMP.dual.(m.ext[:constraints][:con13])
            λ_bc_14 = JuMP.dual.(m.ext[:constraints][:con14])
            λ_bc_15 = JuMP.dual.(m.ext[:constraints][:con15])
        end
    elseif read_data == "5-node_HM" || read_data == "5-node_OBZ"
        λ_bc_7 = JuMP.dual.(m.ext[:constraints][:con7]) 
    end    
end


TOT_opt = zeros(length(T))


# Optimal values for redispatch outcome 
global redispatch_needed = zeros(length(a.ext[:sets][:T])) 
global overload = zeros(length(a.ext[:sets][:T]),length(L))
global overload_DC = zeros(length(a.ext[:sets][:T]),length(L_DC))
for l in 1:length(a.ext[:sets][:L])
    for t in 1:length(a.ext[:sets][:T])

    if abs(a.ext[:parameters][:F_DA][l,t]) > a.ext[:parameters][:TC][l] 
        global redispatch_needed[t] = 1
        global overload[t,l] = 1
        # @info "There is congestion on an AC-line l.", l
    end 

    end
end 
for l_dc in 1:length(a.ext[:sets][:L_DC])
    for t in 1:length(a.ext[:sets][:T])
        
    if abs(a.ext[:parameters][:F_DC_DA][l_dc,t]) > a.ext[:parameters][:TC_DC][l_dc] 
        global redispatch_needed[t] = 1
        global overload_DC[t,l_dc] = 1
        # @info "There is congestion on a DC-line l_dc.", l_dc
     end 

    end
end 

# congestion_DC = zeros(length(a.ext[:sets][:T]), length(L_DC))
# # Loop over all DC lines and all hours
# for l_dc in 1:length(a.ext[:sets][:L_DC])
#     for t in 1:length(a.ext[:sets][:T])
#         # Check if F_DC_DA equals TC_DC
#         if abs(a.ext[:parameters][:F_DC_DA][l_dc, t]) == a.ext[:parameters][:TC_DC][l_dc]
#             # If they are equal, set congestion_DC to 1 for that hour
#             congestion_DC[t, l_dc] = 1
#         end
#     end
# end

congestion_DC = zeros(length(a.ext[:sets][:T]), length(L_DC))

# Loop over all hours
for t in 1:length(a.ext[:sets][:T])
    for l_dc in 1:length(a.ext[:sets][:L_DC])
        if abs(a.ext[:parameters][:F_DC_DA][l_dc, t]) == a.ext[:parameters][:TC_DC][l_dc]
            # If they are equal, set congestion_DC to 1 for that hour and DC line
            congestion_DC[t, l_dc] = 1
        end
    end
end


df_congestion_DC = DataFrame(Transpose(Matrix(congestion_DC)), :auto)
CSV.write(joinpath(base_folder_results, "Congestion_DC.csv"), df_congestion_DC)

if redispatch == "yes"
    UP_opt = value.(b.ext[:variables][:UP])
    DOWN_opt = value.(b.ext[:variables][:DOWN])
    F_beforeRD_opt = value.(b.ext[:expressions][:F_beforeRD])
    F_afterRD_opt = value.(b.ext[:expressions][:F_afterRD])
    RD_per_zone_opt = value.(b.ext[:expressions][:RDC_per_zone])
    curt_delta_opt = value.(b.ext[:variables][:curt_delta])
    F_DC_delta_opt = value.(b.ext[:variables][:F_DC_delta])
    if hydrogen == "yes"
        UP_H2_opt = value.(b.ext[:variables][:UP_H2])
        DOWN_H2_opt = value.(b.ext[:variables][:DOWN_H2])
    end

    for t in 1:length(a.ext[:sets][:T])
        if redispatch_needed[t] == 1
        else
            UP_opt[:,t] = zeros(length(G))
            DOWN_opt[:,t] = zeros(length(G))
            F_afterRD_opt[:,t] = F_beforeRD_opt[:,t]
            curt_delta_opt[:,t] = zeros(length(N))
            F_DC_delta_opt[:,t] = zeros(length(L_DC))
            RDC_opt[t] = 0.0
        end 
    
    TOT_opt[t] = GC_opt[t] + RDC_opt[t]
    end
end

#Congestion rent function
function calculate_congestion_rents(F_DC_DA, df_node, df_DC, λ_opt)
    bus_to_zone = Dict(df_node[i, :BusID] => df_node[i, :Zone] for i in 1:size(df_node, 1))
    rents_data = DataFrame(Line = Int[], Hour = Int[], FromZone = Int[], ToZone = Int[], Flow = Float64[], PriceDifference = Float64[], Rent = Float64[])

    # Identify cross-border DC lines and their zone mappings
    for l_dc in 1:size(df_DC, 1)
        from_bus = df_DC[l_dc, :FromBus]
        to_bus = df_DC[l_dc, :ToBus]
        from_zone = bus_to_zone[from_bus]
        to_zone = bus_to_zone[to_bus]

        if from_zone != to_zone
            for t in 1:size(F_DC_DA, 2)  # Assuming columns in F_DC_DA represent time periods
                flow = F_DC_DA[l_dc, t]
                price_diff = λ_opt[to_zone, t] - λ_opt[from_zone, t]  # Assuming price zones correspond to rows in λ_opt
                rent = flow * price_diff
                push!(rents_data, (l_dc, t, from_zone, to_zone, flow, price_diff, rent))
            end
        end
    end
    return rents_data
end
total_rents_data = calculate_congestion_rents(F_DC_DA, df_node, df_DC, λ_opt)


if cap_calc == "yes"
    # Calculate the congestion costs
    p_domain = a.ext[:parameters][:p_domain] = value.(x.ext[:variables][:p])
    p_FB_domain = a.ext[:parameters][:p_FB_domain] = value.(x.ext[:variables][:p_FB])
    F_FBMC_domain = a.ext[:parameters][:F_FBMC_domain] = value.(x.ext[:variables][:F_FBMC])
    F_DC_domain = a.ext[:parameters][:F_DC_domain] = value.(x.ext[:variables][:F_DC])
end    



# Print the optimal values for decision variables
# println("Optimal values for decision variables:")
# println("v_opt = ", v_opt)
# println("p_opt = ", p_opt)
# println("curt_opt = ", curt_opt)
# println("F_AC_opt = ", F_AC_opt)
# if select_model == "Exact_HC"
#     println("F_bar_opt = ", F_bar_opt)
# end    

# Print the optimal values for dual variables
# println("Optimal values for dual variables:")
# println("λ_opt = ", λ_opt)
# println("λ_opt2 = ", λ_opt2)
# println("λ_opt4 = ", λ_opt4)

# Print the optimal values for computed outcomes
# println("Optimal values for computed outcomes:")
# println("GC_opt = ", GC_opt)



#WRITING THE FILES
CSV.write(joinpath(base_folder_results, "Congestion_Rents.csv"), total_rents_data)

dfNPTDF = DataFrame(NPTDF_opt,:auto)
CSV.write(joinpath(base_folder_results, "NPTDF.csv"), dfNPTDF)

if redispatch == "yes"
    df1 = DataFrame(cost_generation = GC_opt, cost_redispatch = RDC_opt, total_cost = TOT_opt)
    CSV.write(joinpath(base_folder_results, "redispatch/Cost_total.csv"), df1)

    dfa = DataFrame(RDC_perzone_opt.data,:auto)
    CSV.write(joinpath(base_folder_results, "redispatch/Cost_redispatch_per_zone.csv"), dfa)

    dfb = DataFrame(UP_opt.data,:auto)
    CSV.write(joinpath(base_folder_results, "redispatch/Redispatch_Upward.csv"), dfb)

    dfc = DataFrame(DOWN_opt.data,:auto)
    CSV.write(joinpath(base_folder_results, "redispatch/Redispatch_Downward.csv"), dfc)

    dfd = DataFrame(F_beforeRD_opt.data,:auto)
    CSV.write(joinpath(base_folder_results, "redispatch/FlowBeforeRD.csv"), dfd)

    dfe = DataFrame(F_afterRD_opt.data,:auto)
    CSV.write(joinpath(base_folder_results, "redispatch/FlowAfterRD.csv"), dfe)
        
    # df13 = DataFrame(b.ext[:variables][:curt_delta].data,:auto)
    # CSV.write(joinpath(base_folder_results, "redispatch/Redispatch_curtailment.csv"), df13)
    
    df13b = DataFrame(Matrix(curt_delta_opt),:auto)
    CSV.write(joinpath(base_folder_results, "redispatch/Redispatch_curtailment_opt.csv"), df13b)

    # df14 = DataFrame(b.ext[:variables][:F_DC_delta].data,:auto)
    # CSV.write(joinpath(base_folder_results, "redispatch/Redispatch_DCflow.csv"), df14)

    df14b = DataFrame(Matrix(F_DC_delta_opt),:auto)
    CSV.write(joinpath(base_folder_results, "redispatch/DC_delta_RD.csv"), df14b)

    dfλ_rd_1 = DataFrame(shadpw_price_con1 = λ_rd_1.data)
    CSV.write(joinpath(base_folder_results, "redispatch/shadow_price_con1.csv"), dfλ_rd_1)

    dfλ_rd_2 = DataFrame(λ_rd_2.data,:auto)
    CSV.write(joinpath(base_folder_results, "redispatch/shadow_price_con2.csv"), dfλ_rd_2)

    dfλ_rd_3 = DataFrame(λ_rd_3.data,:auto)
    CSV.write(joinpath(base_folder_results, "redispatch/shadow_price_con3.csv"), dfλ_rd_3)

    dfλ_rd_4 = DataFrame(λ_rd_4.data,:auto)
    CSV.write(joinpath(base_folder_results, "redispatch/shadow_price_con4.csv"), dfλ_rd_4)

    dfλ_rd_5 = DataFrame(λ_rd_5.data,:auto)
    CSV.write(joinpath(base_folder_results, "redispatch/shadow_price_con5.csv"), dfλ_rd_5)

    dfλ_rd_6 = DataFrame(λ_rd_6.data,:auto)
    CSV.write(joinpath(base_folder_results, "redispatch/shadow_price_con6.csv"), dfλ_rd_6)

    dfλ_rd_7 = DataFrame(λ_rd_7.data,:auto)
    CSV.write(joinpath(base_folder_results, "redispatch/shadow_price_con7.csv"), dfλ_rd_7)

    dfλ_rd_8 = DataFrame(λ_rd_8.data,:auto)
    CSV.write(joinpath(base_folder_results, "redispatch/shadow_price_con8.csv"), dfλ_rd_8)

    if read_data == "5-node_HM" || read_data == "5-node_OBZ"
        dfλ_rd_9 = DataFrame(shadpw_price_con9 = λ_rd_9.data)
        CSV.write(joinpath(base_folder_results, "redispatch/shadow_price_con9.csv"), dfλ_rd_9)   

    elseif read_data == "Schonheit_HM" || read_data == "Schonheit_OBZ"  || read_data == "Schonheit_OBZ_adjusted" || read_data == "Reference_Case" || read_data == "Simple_Hybrid"
        dfλ_rd_9 = DataFrame(shadpw_price_con9 = λ_rd_9.data)
        CSV.write(joinpath(base_folder_results, "redispatch/shadow_price_con9.csv"), dfλ_rd_9)

        dfλ_rd_10 = DataFrame(shadpw_price_con10 = λ_rd_10.data)
        CSV.write(joinpath(base_folder_results, "redispatch/shadow_price_con10.csv"), dfλ_rd_10)

        dfλ_rd_11 = DataFrame(shadpw_price_con11 = λ_rd_11.data)
        CSV.write(joinpath(base_folder_results, "redispatch/shadow_price_con11.csv"), dfλ_rd_11)

        dfλ_rd_12 = DataFrame(shadpw_price_con12 = λ_rd_12.data)
        CSV.write(joinpath(base_folder_results, "redispatch/shadow_price_con12.csv"), dfλ_rd_12)

        dfλ_rd_13 = DataFrame(shadpw_price_con13 = λ_rd_13.data)
        CSV.write(joinpath(base_folder_results, "redispatch/shadow_price_con13.csv"), dfλ_rd_13)

        if read_data == "Schonheit_HM" || read_data == "Schonheit_OBZ" || read_data == "Schonheit_OBZ_adjusted" || read_data == "Reference_Case"
            dfλ_rd_14 = DataFrame(shadpw_price_con14 = λ_rd_14.data)
            CSV.write(joinpath(base_folder_results, "redispatch/shadow_price_con14.csv"), dfλ_rd_14)
        end
    end

    else
    df1 = DataFrame(cost_generation = GC_opt)
    CSV.write(joinpath(base_folder_results, "Cost_total.csv"), df1)  
end

df21 = DataFrame(RES_prod,:auto)
res_prod_t = DataFrame(transpose(Matrix(df21)), :auto)
CSV.write(joinpath(base_folder_results, "RES_production.csv"), res_prod_t)

if select_model == "Exact_HC"
    df111 = DataFrame(cost_generation = GC_opt)
    CSV.write(joinpath(base_folder_results, "Cost_total.csv"), df1)
        
    df4b = DataFrame(λ_opt4.data,:auto)
    CSV.write(joinpath(base_folder_results, "Price_Exact.csv"), df4b)
       
    df3 = DataFrame(vbar_opt,:auto)
    CSV.write(joinpath(base_folder_results, "Dispatch_bar.csv"), df3)

    df7 = DataFrame(curtbar_opt,:auto)
    CSV.write(joinpath(base_folder_results, "curtbar_opt.csv"), df7)

    df23 = DataFrame(F_bar_opt,:auto)
    CSV.write(joinpath(base_folder_results, "Flow_AC_bar.csv"), df23)

    df24 = DataFrame(F_DC_DA_bar.data,:auto)
    CSV.write(joinpath(base_folder_results, "Flow_DC_bar.csv"), df24)
        
else
    df8 = DataFrame(p_FB_DA,:auto)
    CSV.write(joinpath(base_folder_results, "FlowBased_Position.csv"), df8)

    df9 = DataFrame(F_FBMC_DA,:auto)
    CSV.write(joinpath(base_folder_results, "Commercial_Flow.csv"), df9)    

    df_bc1 = DataFrame(GC_bc,:auto)
    CSV.write(joinpath(base_folder_results, "base-case/Cost_generation_base_case.csv"), df_bc1)
        
    df_bc2 = DataFrame(v_bc,:auto)
    CSV.write(joinpath(base_folder_results, "base-case/Dispatch_base_case.csv"), df_bc2)

    df_bc3 = DataFrame(curt_bc.data, :auto)
    CSV.write(joinpath(base_folder_results, "base-case/Curtailment_base_case.csv"), df_bc3)

    net_production_bc = res_prod_t .- df_bc3
    df_net_production_bc = joinpath(base_folder_results, "base-case/Net_Production_OWF_base_case.csv")
    CSV.write(df_net_production_bc, net_production_bc)
        
    df_bc4 = DataFrame(p_bc,:auto)
    CSV.write(joinpath(base_folder_results, "base-case/NetPosition_base_case.csv"), df_bc4)

    # for i in 1:size(df_bc4, 1)
    #     filename = joinpath(base_folder_results, "FB_domain/NetPosition_bc_z$(i).csv")
    #     row_df = view(df_bc4, i:i, :)
    #     CSV.write(filename, row_df)
    # end

    df_bc5 = DataFrame(F_bc,:auto)
    CSV.write(joinpath(base_folder_results, "base-case/Flow_AC_base_case.csv"), df_bc5)

    df_bc6 = DataFrame(F_dc_bc,:auto)
    CSV.write(joinpath(base_folder_results, "base-case/Flow_DC_base_case.csv"), df_bc6)

    dfλ_bc_1 = DataFrame(λ_bc_1.data,:auto)
    CSV.write(joinpath(base_folder_results, "base-case/shadow_price_con1.csv"), dfλ_bc_1)

    dfλ_bc_2 = DataFrame(λ_bc_2.data,:auto)
    CSV.write(joinpath(base_folder_results, "base-case/shadow_price_con2.csv"), dfλ_bc_2)

    dfλ_bc_3 = DataFrame(λ_bc_3.data,:auto)
    CSV.write(joinpath(base_folder_results, "base-case/shadow_price_con3.csv"), dfλ_bc_3)

    dfλ_bc_4 = DataFrame(λ_bc_4.data,:auto)
    CSV.write(joinpath(base_folder_results, "base-case/shadow_price_con4.csv"), dfλ_bc_4)

    dfλ_bc_5 = DataFrame(λ_bc_5.data,:auto)
    CSV.write(joinpath(base_folder_results, "base-case/shadow_price_con5.csv"), dfλ_bc_5)

    dfλ_bc_6 = DataFrame(λ_bc_6.data,:auto)
    CSV.write(joinpath(base_folder_results, "base-case/shadow_price_con6.csv"), dfλ_bc_6)

    if read_data == "5-node_HM" || read_data == "5-node_OBZ"
        dfλ_bc_7 = DataFrame(Shadow_Price_Con7 = λ_bc_7.data)
        CSV.write(joinpath(base_folder_results, "base-case/shadow_price_con7.csv"), dfλ_bc_7)
    elseif read_data == "Schonheit_HM" || read_data == "Schonheit_OBZ" || read_data == "Schonheit_OBZ_adjusted" || read_data == "Reference_Case" || read_data == "Simple_Hybrid"
        dfλ_bc_7 = DataFrame(Shadow_Price_Con7 = λ_bc_7.data)
        CSV.write(joinpath(base_folder_results, "base-case/shadow_price_con7.csv"), dfλ_bc_7)

        dfλ_bc_8 = DataFrame(Shadow_Price_Con8 = λ_bc_8.data)
        CSV.write(joinpath(base_folder_results, "base-case/shadow_price_con8.csv"), dfλ_bc_8)

        dfλ_bc_9 = DataFrame(ShadowPrice9 = λ_bc_9.data)
        CSV.write(joinpath(base_folder_results, "base-case/shadow_price_con9.csv"), dfλ_bc_9)

        dfλ_bc_10 = DataFrame(ShadowPrice10 = λ_bc_10.data)
        CSV.write(joinpath(base_folder_results, "base-case/shadow_price_con10.csv"), dfλ_bc_10)

        dfλ_bc_11 = DataFrame(ShadowPrice11 = λ_bc_11.data)
        CSV.write(joinpath(base_folder_results, "base-case/shadow_price_con11.csv"), dfλ_bc_11)

        dfλ_bc_12 = DataFrame(ShadowPrice12 = λ_bc_12.data)
        CSV.write(joinpath(base_folder_results, "base-case/shadow_price_con12.csv"), dfλ_bc_12)

        if read_data == "Schonheit_HM" || read_data == "Schonheit_OBZ" || read_data == "Schonheit_OBZ_adjusted" || read_data == "Reference_Case"
            dfλ_bc_13 = DataFrame(ShadowPrice12 = λ_bc_13.data)
            CSV.write(joinpath(base_folder_results, "base-case/shadow_price_con13.csv"), dfλ_bc_13)

            dfλ_bc_14 = DataFrame(ShadowPrice12 = λ_bc_14.data)
            CSV.write(joinpath(base_folder_results, "base-case/shadow_price_con14.csv"), dfλ_bc_14)

            dfλ_bc_15 = DataFrame(ShadowPrice12 = λ_bc_15.data)
            CSV.write(joinpath(base_folder_results, "base-case/shadow_price_con15.csv"), dfλ_bc_15)
        end
    end

    # dfrp = DataFrame(RAM_pos_opt,:auto)
    # CSV.write(joinpath(base_folder_results, "RAM_pos.csv"), dfrp)

    # for i in 1:size(dfrp, 1)
    #     filename = joinpath(base_folder_results, "FB_domain/RAM_pos_l$(i).csv")
    #     row_df = view(dfrp, i:i, :)
    #     CSV.write(filename, row_df)
    # end

    # dfrn = DataFrame(RAM_neg_opt,:auto)
    # CSV.write(joinpath(base_folder_results, "RAM_neg.csv"), dfrn)

    # for i in 1:size(dfrn, 1)
    #     filename = joinpath(base_folder_results, "FB_domain/RAM_neg_l$(i).csv")
    #     row_df = view(dfrn, i:i, :)
    #     CSV.write(filename, row_df)
    # end

    # dfZPTDF = DataFrame(ZPTDF_opt,:auto)
    # CSV.write(joinpath(base_folder_results, "ZPTDF.csv"), dfZPTDF)

    # for i in 1:size(dfZPTDF, 1)
    #     filename = joinpath(base_folder_results, "FB_domain/ZPTDF_z$(i).csv")
    #     row_df = view(dfZPTDF, i:i, :)
    #     CSV.write(filename, row_df)
    # end

    dfrp_bc = DataFrame(RAM_pos_opt_bc,:auto)
    CSV.write(joinpath(base_folder_results, "base-case/RAM_pos_base_case.csv"), dfrp_bc)
    
    # for i in 1:size(dfrp_bc, 1)
    #     filename = joinpath(base_folder_results, "FB_domain/RAM_pos_bc_l$(i).csv")
    #     row_df = view(dfrp_bc, i:i, :)
    #     CSV.write(filename, row_df)
    # end

    dfrn_bc = DataFrame(RAM_neg_opt_bc,:auto)
    CSV.write(joinpath(base_folder_results, "base-case/RAM_neg_base_case.csv"), dfrn_bc)

    # for i in 1:size(dfrn_bc, 1)
    #     filename = joinpath(base_folder_results, "FB_domain/RAM_neg_bc_l$(i).csv")
    #     row_df = view(dfrn_bc, i:i, :)
    #     CSV.write(filename, row_df)
    # end

    dfZPTDF_bc = DataFrame(ZPTDF_opt_bc,:auto)
    CSV.write(joinpath(base_folder_results, "base-case/ZPTDF_base_case.csv"), dfZPTDF_bc)

    # for i in 1:size(dfZPTDF_bc, 1)
    #     filename = joinpath(base_folder_results, "FB_domain/ZPTDF_bc_z$(i).csv")
    #     row_df = view(dfZPTDF_bc, i:i, :)
    #     CSV.write(filename, row_df)
    # end
end
   
df2 = DataFrame(GC_DA.data,:auto)
CSV.write(joinpath(base_folder_results, "Cost_generation.csv"), df2)
    
df4 = DataFrame(λ_opt.data,:auto)
CSV.write(joinpath(base_folder_results, "Price.csv"), df4)
   
df5 = DataFrame(v_opt.data,:auto)
CSV.write(joinpath(base_folder_results, "Dispatch.csv"), df5)
   
df6 = DataFrame(p_opt.data,:auto)
CSV.write(joinpath(base_folder_results, "Total_NetPosition.csv"), df6) 

# for i in 1:size(df6, 1)
#     filename = joinpath(base_folder_results, "FB_domain/NetPosition_z$(i).csv")
#     row_df = view(df6, i:i, :)
#     CSV.write(filename, row_df)
# end

df11 = DataFrame(curt_DA.data,:auto)
CSV.write(joinpath(base_folder_results, "Curtailment.csv"), df11)
   
net_production = res_prod_t .- df11
df_net_production = joinpath(base_folder_results, "Net_production_OWF.csv")
CSV.write(df_net_production, net_production)

df12 = DataFrame(F_DC_DA.data,:auto)
CSV.write(joinpath(base_folder_results, "Flow_DC.csv"), df12)
  
df15 = DataFrame(overload,:auto)
CSV.write(joinpath(base_folder_results, "Overload.csv"), df15)
   
df16 = DataFrame(overload_DC,:auto)
CSV.write(joinpath(base_folder_results, "Overload_DC.csv"), df16)
    
df17 = DataFrame(i,:auto)
CSV.write(joinpath(base_folder_results, "i.csv"), df17)
    
df18 = DataFrame(i_DC,:auto)
CSV.write(joinpath(base_folder_results, "i_DC.csv"), df18)
    

    
df22 = DataFrame(F_AC_opt,:auto)
CSV.write(joinpath(base_folder_results, "Flow_AC.csv"), df22) 

if hydrogen == "yes"
    dfh2 = DataFrame(e_DA,:auto)
    CSV.write(joinpath(base_folder_results, "Hydrogen_dispatch.csv"), dfh2)

    dfh2_bc = DataFrame(e_bc,:auto)
    CSV.write(joinpath(base_folder_results, "base-case/Hydrogen_dispatch_base_case.csv"), dfh2_bc)

    dfh2_UP = DataFrame(UP_H2_opt.data,:auto)
    CSV.write(joinpath(base_folder_results, "redispatch/Hydrogen_Upward.csv"), dfh2_UP)

    dfh2_DOWN = DataFrame(DOWN_H2_opt.data,:auto)
    CSV.write(joinpath(base_folder_results, "redispatch/Hydrogen_Downward.csv"), dfh2_DOWN)

    dfλ6 = DataFrame(shadow_price_con6 = λ_opt6.data)
    CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con6.csv"), dfλ6)

    dfh2λ_rd_h2_1 = DataFrame(λ_rd_h2_1.data, :auto)
    CSV.write(joinpath(base_folder_results, "redispatch/shadow_price_conh2_1.csv"), dfh2λ_rd_h2_1)

    dfh2λ_rd_h2_2 = DataFrame(λ_rd_h2_2.data, :auto)
    CSV.write(joinpath(base_folder_results, "redispatch/shadow_price_conh2_2.csv"), dfh2λ_rd_h2_2)
end
    
# #Code below is to analyse the shadow prices of the dual variables
dfλ = DataFrame(λ_opt.data,:auto)
CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con1.csv"), dfλ)
   
dfλ2 = DataFrame(λ_opt2.data,:auto)
CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con2.csv"), dfλ2)

dfλ3 = DataFrame(λ_opt3.data,:auto)
CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con3.csv"), dfλ3)

dfλ4 = DataFrame(λ_opt4.data,:auto)
CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con4.csv"), dfλ4)

dfλ5 = DataFrame(λ_opt5.data,:auto)
CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con5.csv"), dfλ5)



if read_data == "5-node_HM" || read_data == "5-node_OBZ"
    dfλ6 = DataFrame(shadow_price_con6 = λ_opt6.data)
    CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con6.csv"), dfλ6)
   
    dfλ7 = DataFrame(Shadow_Price_Con7 = λ_opt7.data)
    CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con7.csv"), dfλ7)

    dfλ8 = DataFrame(Shadow_Price_Con8 = λ_opt8.data)
    CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con8.csv"), dfλ8)

elseif read_data == "Schonheit_HM" || read_data == "Schonheit_OBZ" || read_data == "Schonheit_OBZ_adjusted" || read_data == "Reference_Case" || read_data == "Simple_Hybrid"
    dfλ6 = DataFrame(shadow_price_con6 = λ_opt6.data)
    CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con6.csv"), dfλ6)
    
    dfλ7 = DataFrame(Shadow_Price_Con7 = λ_opt7.data)
    CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con7.csv"), dfλ7)
    
    dfλ8 = DataFrame(Shadow_Price_Con8 = λ_opt8.data)
    CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con8.csv"), dfλ8)

    dfλ9 = DataFrame(ShadowPrice9 = λ_opt9.data)
    CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con9.csv"), dfλ9)
    
    dfλ10 = DataFrame(ShadowPrice10 = λ_opt10.data)
    CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con10.csv"), dfλ10)

    dfλ11 = DataFrame(ShadowPrice11 = λ_opt11.data)
    CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con11.csv"), dfλ11)

    dfλ12 = DataFrame(ShadowPrice12 = λ_opt12.data)
    CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con12.csv"), dfλ12)

    dfλ13 = DataFrame(ShadowPrice13 = λ_opt13.data)
    CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con13.csv"), dfλ13)

    dfλ14 = DataFrame(ShadowPrice14 = λ_opt14.data)
    CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con14.csv"), dfλ14)

    dfλ15 = DataFrame(ShadowPrice15 = λ_opt15.data)
    CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con15.csv"), dfλ15)

    dfλ16 = DataFrame(ShadowPrice16 = λ_opt16.data)
    CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con16.csv"), dfλ16)
    
    if read_data == "Schonheit_HM" || read_data == "Schonheit_OBZ" || read_data == "Schonheit_OBZ_adjusted" || read_data == "Reference_Case"

        dfλ17 = DataFrame(ShadowPrice17 = λ_opt17.data)
        CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con17.csv"), dfλ17)

        dfλ18 = DataFrame(ShadowPrice18 = λ_opt18.data)
        CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con18.csv"), dfλ18)

        dfλ19 = DataFrame(ShadowPrice19 = λ_opt19.data)
        CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con19.csv"), dfλ19)

        dfλ20 = DataFrame(ShadowPrice20 = λ_opt20.data)
        CSV.write(joinpath(base_folder_results, "shadow_prices/shadow_price_con20.csv"), dfλ20)
        
    end
end

if cap_calc == "yes"
    dfFBd2 = DataFrame(p_domain.data,:auto)
    CSV.write(joinpath(base_folder_results, "FB_domain/MaxOBZ_NPs.csv"), dfFBd2)

    dfFBd3 = DataFrame(p_FB_domain.data,:auto)
    CSV.write(joinpath(base_folder_results, "FB_domain/Flow_based_positions.csv"), dfFBd3)

    dfFB4 = DataFrame(F_FBMC_domain.data,:auto)
    CSV.write(joinpath(base_folder_results, "FB_domain/max_possible_AC_flows.csv"), dfFB4)

    dfFB5 = DataFrame(F_DC_domain.data,:auto)
    CSV.write(joinpath(base_folder_results, "FB_domain/max_possible_DC_flows.csv"), dfFB5)
end

flow_FBMC_DA = a.ext[:parameters][:Flow_FBMC]
df_flow_FBMC_DA = DataFrame(flow_FBMC_DA.data,:auto)
CSV.write(joinpath(base_folder_results, "flow_FBMC_DA.csv"), df_flow_FBMC_DA)


