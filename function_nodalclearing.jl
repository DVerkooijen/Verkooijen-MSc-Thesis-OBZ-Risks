m = Model(Gurobi.Optimizer)


function market_clearing!(a::Model)
    a.ext[:variables] = Dict()      # To store variables, expressions, constraints and objective in Dictionary
    a.ext[:expressions] = Dict()
    a.ext[:constraints] = Dict()
    a.ext[:objective] = Dict()

    # Extract sets
    N = a.ext[:sets][:N]        #Extracted from 'external' data file
    G = a.ext[:sets][:G]
    Z = a.ext[:sets][:Z]
    L = a.ext[:sets][:L]
    L_DC = a.ext[:sets][:L_DC]
    T = a.ext[:sets][:T]

    # Extract parameters
    GENCAP = a.ext[:parameters][:GENCAP]    #Extracted from 'external' data file
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
    v = a.ext[:variables][:v] = @variable(a, [g in G, t in T], lower_bound = 0, upper_bound = 1, base_name = "dispatch")
    curt = a.ext[:variables][:curt] = @variable(a, [n in N, t in T], lower_bound = 0, base_name = "curtailment")
    p = a.ext[:variables][:p] = @variable(a, [n in N, t in T], base_name = "position")
    F = a.ext[:variables][:F] = @variable(a, [l in L, t in T], base_name = "flows")
    F_dc = a.ext[:variables][:F_dc] = @variable(a, [l_dc in L_DC, t in T], base_name = "DC flows")

    # Create expression Congestion Costs per time
    GC_per_time_and_zone = a.ext[:expressions][:GC_per_zone] = @expression(a, [z in Z, t in T], sum(sum(MC[g]*GENCAP[g]*v[g,t] for g in g_in_n[n]) + CC*curt[n,t] for n in n_in_z[z]) )
    
    # Objective
    GC = a.ext[:objective][:GC] = @objective(a, Min, sum(sum(MC[g]*GENCAP[g]*v[g,t] for g in g_in_n[n]) + CC*curt[n,t] for t in T for n in N))

    # Constraints 
    
    #Power Balance for each generator n at each time t
    a.ext[:constraints][:con1] = @constraint(a, con1[n in N, t in T], sum(GENCAP[g] * v[g,t] for g in g_in_n[n]) + -curt[n,t] - DEM[t,findfirst(N .== n)] - p[n,t] == 0 )
    #power balance of each load n at each time t
    a.ext[:constraints][:con2] = @constraint(a, con2[n in N, t in T], sum(F[l,t]*inc[l,findfirst(N .== n)] for l in L) - p[n,t] == 0 )
    #power flow on each line l at each time t (Based on NPTDF)   
    a.ext[:constraints][:con3] = @constraint(a, con3[l in L, t in T], F[l,t] == sum(NPTDF[findfirst(N .== n),l]*(sum(GENCAP[g]*(v[g,t]) for g in g_in_n[n]) -curt[n,t] -DEM[t,findfirst(N .== n)]) for n in N ))
    #upper and lower bounds of power flows for each tranmission line
    a.ext[:constraints][:con4] = @constraint(a, con4[l in L, t in T], -TC[l] <= F[l,t] <= TC[l] )
    #Curtailed amount of each generator n at each time t is less or equal to available Ren energy
    a.ext[:constraints][:con6] = @constraint(a, con6[n in N, t in T], curt[n,t] <= sum(GENCAP[g] * v[g, t] for g in g_in_n[1]))

    return a
end

# Build your model
market_clearing!(a)

optimize!(a)

#IF the model gives an infeasible solution
if termination_status(a) == MOI.INFEASIBLE
    println("Model is infeasible. Computing IIS.")
    
    # Qualify the use of `backend` to resolve the namespace conflict
    optimizer_model = JuMP.backend(a).optimizer.model
    
    # Access the Gurobi model directly
    gurobi_model = optimizer_model.inner

    # Use the Gurobi function to compute the IIS
    Gurobi.GRBcomputeIIS(gurobi_model)

    # Optionally, save the IIS information to a file
    Gurobi.GRBwrite(gurobi_model, "model.ilp")
end

#Extract the values of the variables and expressions after optimization
GC_DA = a.ext[:parameters][:GC_DA] = value.(a.ext[:expressions][:GC_per_zone])
v_DA = a.ext[:parameters][:v_DA] = value.(a.ext[:variables][:v])
curt_DA = a.ext[:parameters][:curt_DA] = value.(a.ext[:variables][:curt])
p_DA = a.ext[:parameters][:p_DA] = value.(a.ext[:variables][:p])
F_DA = a.ext[:parameters][:F_DA] = value.(a.ext[:variables][:F])
F_DC_DA = a.ext[:parameters][:F_DC_DA] = value.(a.ext[:variables][:F_dc])
RES_prod = a.ext[:parameters][:RES_prod] = value.(a.ext[:parameters][:RES])

println("Optimal value of the objective function: ", objective_value(a))