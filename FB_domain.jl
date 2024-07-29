using JuMP
using Gurobi
using Printf
using CSV
using DataFrames
using Plots

# Parameters
base_folder_results = a.ext[:loops][:base_folder_results]
hour = 343          #Change desried hour to take snapshot
zone_x_axis = 1     #change zones to take as x and y axis
zone_y_axis = 2

max_OBZ = "no"      #yes if OBZ optimisation is to be inlcuded (no = regular FB domain)
snapshot ="first"   #Change for saving purposes of the png (if we have 2 graphs desribing 1 snapshot, there is first to sixt graph inclduded)
D1_clearing = "no"  #yes if D-1 clearing point is to be included

# Load parameters
RAM_pos_bc = a.ext[:parameters][:RAM_pos_bc]
RAM_neg_bc = a.ext[:parameters][:RAM_neg_bc]
RAM_pos_max = a.ext[:parameters][:RAM_pos_max]
RAM_neg_max = a.ext[:parameters][:RAM_neg_max]

ZPTDF = a.ext[:parameters][:ZPTDF_bc]
cross_border_lines = CSV.read(joinpath(base_folder_results, "FB_domain/cross_border_lines.csv"), DataFrame)
ntc_values = CSV.read(joinpath(base_folder_results, "FB_domain/ntc_values.csv"), DataFrame)

function compute_line_parameters(ram, ptdf_values)
    intercept = ram
    return ptdf_values, intercept
end

struct Formula
    line::Int
    type::String
    a::Float64
    b::Float64
    c::Float64
end

#this is with the RAM*ZPTDF
function generate_linear_formulas(RAM_pos, RAM_neg, ZPTDF, ntc_values, cross_border_lines, zone_x_axis, zone_y_axis, hours)
    formulas = Formula[]
    
    for hour in hours
        # Calculate the sum of RAM * |ZPTDF| for all AC lines connected to the zone
        ac_capacity_sum_x = 0.0
        ac_capacity_sum_y = 0.0
        for ac_row in eachrow(cross_border_lines)
            if ac_row.LineType == "AC"
                branch_id = ac_row.BranchID
                if ac_row.FromZone == zone_x_axis || ac_row.ToZone == zone_x_axis
                    ram_value = RAM_pos[branch_id, hour]  # Adjust to use RAM_neg if needed for negative scenarios
                    ac_capacity_sum_x += ram_value
                end
                if ac_row.FromZone == zone_y_axis || ac_row.ToZone == zone_y_axis
                    ram_value = RAM_pos[branch_id, hour]  # Adjust to use RAM_neg if needed for negative scenarios
                    ac_capacity_sum_y += ram_value
                end
            end
        end

        # Calculate transmission constraints based on PTDF and RAM for existing lines
        for line_index in 1:size(ZPTDF, 2)
            ram_pos_value = RAM_pos[line_index, hour]
            ram_neg_value = -RAM_neg[line_index, hour]

            ptdf_values = ZPTDF[:, line_index]

            ptdf_coefficients_pos, intercept_pos = compute_line_parameters(ram_pos_value, ptdf_values)
            ptdf_coefficients_neg, intercept_neg = compute_line_parameters(ram_neg_value, ptdf_values)

            push!(formulas, Formula(line_index, "CNE", ptdf_coefficients_pos[zone_x_axis], ptdf_coefficients_pos[zone_y_axis], intercept_pos))
            push!(formulas, Formula(line_index, "CNE", ptdf_coefficients_neg[zone_x_axis], ptdf_coefficients_neg[zone_y_axis], intercept_neg))

            # println("Formula Pos: ", ptdf_coefficients_pos[zone_x_axis], "*x + ", ptdf_coefficients_pos[zone_y_axis], "*y = ", intercept_pos)
            # println("Formula Neg: ", ptdf_coefficients_neg[zone_x_axis], "*x + ", ptdf_coefficients_neg[zone_y_axis], "*y = ", intercept_neg)
        end

       # Adjust constraints for DC lines by incorporating AC line capacity
        for dc_row in eachrow(cross_border_lines)
            if dc_row.LineType == "DC"
                branch_id = dc_row.BranchID
                ntc_row = findfirst(==(branch_id), ntc_values.BranchID)
                ntc_value = ntc_row !== nothing ? ntc_values[ntc_row, :NTC] : 0.0

                if dc_row.FromZone == zone_x_axis || dc_row.ToZone == zone_x_axis
                    # Add constraints for both positive and negative directions including AC capacity
                    push!(formulas, Formula(branch_id, "NTC", 1.0, 0.0, ntc_value + ac_capacity_sum_x))
                    push!(formulas, Formula(branch_id, "NTC", 1.0, 0.0, -(ntc_value + ac_capacity_sum_x)))
                end
                if dc_row.FromZone == zone_y_axis || dc_row.ToZone == zone_y_axis
                    # Add constraints for both positive and negative directions including AC capacity
                    push!(formulas, Formula(branch_id, "NTC", 0.0, 1.0, ntc_value + ac_capacity_sum_y))
                    push!(formulas, Formula(branch_id, "NTC", 0.0, 1.0, -(ntc_value + ac_capacity_sum_y)))
                end
            end
        end

        CSV.write(joinpath(base_folder_results, "FB_domain/formulas.csv"), formulas)
    end
    

    return formulas
end

# Generate formulas for a specific range of hours (example for hour 2)
formulas_bc = generate_linear_formulas(RAM_pos_bc, RAM_neg_bc, ZPTDF, ntc_values, cross_border_lines, zone_x_axis, zone_y_axis, [hour])
formulas_max = generate_linear_formulas(RAM_pos_max, RAM_neg_max, ZPTDF, ntc_values, cross_border_lines, zone_x_axis, zone_y_axis, [hour])

# Function to plot lines using Plots
function plot_lines_plt!(plt, formulas, color, label)
    x_range = -2000:2000  # Adjust this range as needed
    first_line = true
    for formula in formulas
        a, b, c = formula.a, formula.b, formula.c
        # println("Plotting formula: ", a, "*x + ", b, "*y = ", c)
        if b != 0
            y_values = (c .- a .* x_range) ./ b
            mask = (y_values .>= -5000) .& (y_values .<= 5000)  # Adjust y limit
            if first_line
                plot!(plt, x_range[mask], y_values[mask], color = color, linestyle = :solid, label = label)
                first_line = false
            else
                plot!(plt, x_range[mask], y_values[mask], color = color, linestyle = :solid, label = "")
            end
        elseif a != 0
            x_value = c / a
            if x_value >= -5000 && x_value <= 5000  # Adjust x limit
                if first_line
                    plot!(plt, [x_value, x_value], [-5000, 5000], color = color, linestyle = :solid, label = label)
                    first_line = false
                else
                    plot!(plt, [x_value, x_value], [-5000, 5000], color = color, linestyle = :solid, label = "")
                end
            end
        end
    end
end

# Plotting with Plots
function plot_fb_domain_plt(formulas_bc, formulas_max, xlims, ylims)
    plt = plot(xlims = xlims, ylims = ylims, title = "Flow-Based Domain", 
               xlabel = "Net Position Zone $zone_x_axis [MW]", ylabel = "Net Position Zone $zone_y_axis [MW]",
               grid = :on, xticks = -2000:500:2000, yticks = -2000:500:2000,  # Adjust the tick spacing as needed
               gridalpha = 0.5, gridcolor = :black, gridlinestyle = :dash, 
               legend = :topright)  # Set grid line color and style
    
    plot!(plt, [0, 0], ylims, color = :black, linestyle = :solid, label = "")  # Line at x=0
    plot!(plt, xlims, [0, 0], color = :black, linestyle = :solid, label = "")  # Line at y=0
    
    plot_lines_plt!(plt, formulas_bc, :red, "CNEs")
    if max_OBZ == "yes"
        plot_lines_plt!(plt, formulas_max, :purple, "CNEs'")
    end

    if D1_clearing == "yes"
        total_net_positions = CSV.read(joinpath(base_folder_results, "Total_NetPosition.csv"), DataFrame)
        x_value = total_net_positions[zone_x_axis, hour]
        y_value = total_net_positions[zone_y_axis, hour]
        plot!(plt, [x_value], [y_value], seriestype = :scatter, markercolor = :blue, markershape = :diamond, label = "D-1")
    end


    # Save the plot
    mkpath(joinpath(base_folder_results, "FB_domain/graphs"))

    if snapshot == "first"
        if max_OBZ == "yes"
            savefig(plt, joinpath(base_folder_results, "FB_domain/graphs", "FB_Domain_mxOBZ_first_h$(hour).png"))
        elseif max_OBZ == "no"
            savefig(plt, joinpath(base_folder_results, "FB_domain/graphs", "FB_Domain_first_h$(hour).png"))
        end
    elseif snapshot == "second"
        if max_OBZ == "yes"
            savefig(plt, joinpath(base_folder_results, "FB_domain/graphs", "FB_Domain_mxOBZ_second_h$(hour).png"))
        elseif max_OBZ == "no"
            savefig(plt, joinpath(base_folder_results, "FB_domain/graphs", "FB_Domain_second_h$(hour).png"))
        end
    elseif snapshot == "third"
        if max_OBZ == "yes"
            savefig(plt, joinpath(base_folder_results, "FB_domain/graphs", "FB_Domain_mxOBZ_third_h$(hour).png"))
        elseif max_OBZ == "no"
            savefig(plt, joinpath(base_folder_results, "FB_domain/graphs", "FB_Domain_third_h$(hour).png"))
        end
    elseif snapshot == "fourth"
        if max_OBZ == "yes"
            savefig(plt, joinpath(base_folder_results, "FB_domain/graphs", "FB_Domain_mxOBZ_fourth_h$(hour).png"))
        elseif max_OBZ == "no"
            savefig(plt, joinpath(base_folder_results, "FB_domain/graphs", "FB_Domain_fourth_h$(hour).png"))
        end
    elseif snapshot == "fifth"
        if max_OBZ == "yes"
            savefig(plt, joinpath(base_folder_results, "FB_domain/graphs", "FB_Domain_mxOBZ_fifth_h$(hour).png"))
        elseif max_OBZ == "no"
            savefig(plt, joinpath(base_folder_results, "FB_domain/graphs", "FB_Domain_fifth_h$(hour).png"))
        end
    elseif snapshot == "sixth"
        if max_OBZ == "yes"
            savefig(plt, joinpath(base_folder_results, "FB_domain/graphs", "FB_Domain_mxOBZ_sixth_h$(hour).png"))
        elseif max_OBZ == "no"
            savefig(plt, joinpath(base_folder_results, "FB_domain/graphs", "FB_Domain_sixth_h$(hour).png"))
        end
    end
    plt
end

# Define the axis limits
xlims = (-2000, 2000)
ylims = (-2000, 2000)

# Plot the formulas with axis limits
plot_fb_domain_plt(formulas_bc, formulas_max, xlims, ylims)
