# MSc Thesis "Unlocking Offshore Hybrid Projects: Modelling Price and Volume Risks and their Mitigation Measures of Hybrid Offshore Wind Projects in Offshore Bidding Zones under Flow-Based Market Coupling"


## Background
This model was developped for the MSc Thesis "Unlocking Offshore Hybrid Projects: Modelling Price and Volume Risks and their Mitigation Measures of Hybrid Offshore Wind Projects in Offshore Bidding Zones under Flow-Based Market Coupling" by Daan Verkooijen in fulfillment of the the Master of Science degree of Complex Systems Engineering and Management at the Delft University of Technology (Faculty of Technology, Policy and Management). 
The mathemathical model is taken from Kenis et al. (2023), available at https://github.com/kbruninx/OBZvsHM, and extended fit for purpose of the thesis' objectives (See full explanation in the MSc thesis, published at the TU Delft Repository: https://repository.tudelft.nl/).


## How to implement the model for own application
The main input files are in the folder grid_topology. The files describing the electricity system can be adjusted for more realistic grid representations:
    - node_info.csv --> coordinates of nodes in system + attributed bidding zone
    - line_info.csv --> AC lines connected between nodes with max thermal capacities
    - DC_info.csv --> DC lines connected between nodes with max thermal capacities + From/ToDirection to indicate cross-zonal interconnectors
    - incidence.csv --> incidence matrix of all AC lines
    - incidence_dv.csv --> incidence matrix of all DC lines
    - susceptance.csv --> susceptance matrix of all AC lines
    - load_info.csv --> nodal (inelastic) hourly load values
    - offshore.csv, onshore.csv and pv.csv --> hourly capacity factors for renewable generators per zone

Go to the main.jl file and adjust the loops according to the made changes. Important is to adjust the read_data loop based on the (created) input folder for the data. Additionally, a Results folder should be created where the results per simulation are stored to. Default is Reference_Case (Triangular OBZ setup) and Simple_Hybrid (Dual OBZ setup) for the representative grid. 
If a new grid topology is used, make sure to check all files included in main.jl for the required adjustments in the model, indicated with comments. For computational time purposes, it is advised to simulate per month of data and copy and paste per simulation the results from the folder where results are written to, to a desired folder to store the results. The results can then later be aggregated using aggregate_months.py.
Gurobi Optimization Licence is required for the simulation (free for academics https://www.gurobi.com/academia/academic-program-and-licenses/)

After running the main.jl file, monthly results can be aggregated (create folder based on Results/___Example_Folder) whereafter the risk metrics can be simulated (analyse_risk_metrics.py) to determine the frequency and severity of the price and volume risk. Additionaly, FB_domain.jl provides a script for visualisation of the FB domain per hour and create_graphs.py provides a script for price duration curves and volume risk duration curves. 



### References
M. Kenis, E. Delarue, K. Bruninx and F. Dominguez, "Off-shore Bidding Zones Under Flow-Based Market Coupling," 2023 IEEE Belgrade PowerTech, Belgrade, Serbia, 2023, pp. 1-6, doi: 10.1109/PowerTech55446.2023.10202755.

### Contact
LinkedIn:   https://www.linkedin.com/in/daan-verkooijen-9b9679189/
Email:      dverkooijen.zakelijk@hotmail.com