using Pkg
Pkg.activate(@__DIR__) # @__DIR__ = directory this script is in
# Pkg.instantiate() # If a Manifest.toml file exist in the current project, download all the packages declared in that manifest. Else, resolve a set of feasible packages from the Project.toml files and install them.
# Pkg.resolve()  #Resolve and update packages based on Project.toml and Manifest.toml

using CSV
using XLSX
using DataFrames
using YAML
using Printf
using JuMP
using Gurobi
using Cbc
using DelimitedFiles
using Profile

cd("C:/Users/DVerk/Verkooijen-MSc-Thesis-OBZ-Risks/")           #If you've downloaded the repository locally, change this the the exact path of the folder

a = Model(Gurobi.Optimizer)

#Adjust loops to preferences.
a.ext[:loops] = Dict()
read_data = a.ext[:loops][:read_data] = "Reference_Case"        #Used in thesis: Reference_Case for triangle, Simple_Hybrid otherwise (Or test_networks: Schonheit_OBZ, 5-node_OBZ, etc..). If entirely new grid topology is used, the optimization models must be adjusted to to account for these changes (e.g. see comments in function_zonalclearing_GSK_AHC.jl how to do this)
select_model = a.ext[:loops][:select_model] = "GSK_AHC"         #"GSK_AHC" , "Exact_HC", "GSK_SHC"  (See Kenis et al. (2023) for more info on these representations of the DC lines in the FB algoritm. Only GSK_AHC has been adapted for this thesis)
redispatch = a.ext[:loops][:redispatch] = "yes"                 #"yes", "no"       
cap_calc = a.ext[:loops][:cap_calc] = "yes"                     #"yes", "no"    (This is the optimzation problem to distinghuis between the curtailed volume by the capacity calculation and by the capacity allocation)
hydrogen = a.ext[:loops][:hydrogen] = "no"                      #"yes", "no"


#Add here the base folders to store the results for each case to simulate
if read_data == "5-node_HM"
    base_folder_results = a.ext[:loops][:base_folder_results] = "Results/test_models/5-node_HM"
elseif read_data == "5-node_OBZ"
    base_folder_results = a.ext[:loops][:base_folder_results] = "Results/test_models/5-node_OBZ"
elseif read_data == "Schonheit_HM"
    base_folder_results = a.ext[:loops][:base_folder_results] = "Results/test_models/_Schonheit_HM"
elseif read_data == "Schonheit_OBZ"
    base_folder_results = a.ext[:loops][:base_folder_results] = "Results/test_models/_Schonheit_OBZ"
elseif read_data == "Schonheit_OBZ_adjusted"
    base_folder_results = a.ext[:loops][:base_folder_results] = "Results/test_models/_Schonheit_OBZ_adjusted"
elseif read_data == "Reference_Case"
    base_folder_results = a.ext[:loops][:base_folder_results] = "Results/Reference_Case"
elseif read_data == "Simple_Hybrid"
    base_folder_results = a.ext[:loops][:base_folder_results] = "Results/Simple_Hybrid"
end                


#Add/change in this loop the base folders for input data of new cases/grid topologies to simulate!
if read_data == "5-node_HM"
    base_folder_data = a.ext[:loops][:base_folder_data] = "grid_topology/test_networks/5-node_HM"

elseif read_data == "5-node_OBZ"
    base_folder_data = a.ext[:loops][:base_folder_data] = "grid_topology/test_networks/5-node_OBZ"

elseif read_data == "Schonheit_HM"
    base_folder_data = a.ext[:loops][:base_folder_data] = "grid_topology/test_networks/_Schonheit_HM"

elseif read_data == "Schonheit_OBZ"
    base_folder_data = a.ext[:loops][:base_folder_data] = "grid_topology/test_networks/_Schonheit_OBZ"

elseif read_data == "Schonheit_OBZ_adjusted"
    base_folder_data = a.ext[:loops][:base_folder_data] = "grid_topology/test_networks/_Schonheit_OBZ_adjusted"

elseif read_data == "Reference_Case"
    base_folder_data = a.ext[:loops][:base_folder_data] = "grid_topology/Reference_Case"

elseif read_data == "Simple_Hybrid"
    base_folder_data = a.ext[:loops][:base_folder_data] = "grid_topology/Simple_Hybrid"
end 


#---- Here starts the actual simulation -----

#Make sure to check the assumptions in the data file! (MinRAM, α (=cost mark-up redispatch) and NTC values)
include("data.jl")

if hydrogen == "yes"
    include("function_zonalclearing_GSK_AHC_H2.jl")     #if new grid topology is used (e.g. new read_data loop parameter), make sure to adjust the loops in this file accordingly (follow comments function_zonalclearing_GSK_AHC.jl)
    include("redispatch_H2.jl")                         #if new grid topology is used (e.g. new read_data loop parameter), make sure to adjust the loops in this file accordingly (follow comments redispatch.jl)

else    
    include("function_zonalclearing_GSK_AHC.jl")        #if new grid topology is used (e.g. new read_data loop parameter), make sure to adjust the loops in this file accordingly (follow comments)
    include("redispatch.jl")                            #if new grid topology is used (e.g. new read_data loop parameter), make sure to adjust the loops in this file accordingly (follow comments)
end

# if select_model == "Exact_HC"
#     include("function_zonalclearing_exact_HC.jl")
# end


if cap_calc == "yes"
    include("max_net_position_OBZ.jl")                 #if new grid topology is used (e.g. new read_data loop parameter), make sure to adjust the loops in this file accordingly (follow comments)
end

include("process_results.jl")                           #if new grid topology is used (e.g. new read_data loop parameter), make sure to adjust the loops in this file accordingly (primairy changes are the shadow prices outputs)




##These optimization functions have not been included in this thesis, but are available in the repository:
    # include("function_zonalclearing_exact_FBonly.jl")

    # include("function_zonalclearing_GSK_FBonly.jl")

    # include("function_nodalclearing.jl")

    # include("function_zonalclearing_GSK_SHC.jl")