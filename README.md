# MSc Thesis "Unlocking Offshore Hybrid Projects: Modelling Price and Volume Risks and their Mitigation Measures for Offshore Hybrid Projects in Offshore Bidding Zones under Flow-Based Market Coupling"

## Background
This model was developed for the MSc Thesis "Unlocking Offshore Hybrid Projects: Modelling Price and Volume Risks and their Mitigation Measures for Offshore Hybrid Projects in Offshore Bidding Zones under Flow-Based Market Coupling" by Daan Verkooijen, in fulfillment of the Master of Science degree in Complex Systems Engineering and Management at Delft University of Technology (Faculty of Technology, Policy and Management).

The mathematical model is adapted from Kenis et al. (2023), available at [GitHub](https://github.com/kbruninx/OBZvsHM), and extended for the thesis' objectives (see the full explanation in the MSc thesis, published at the [TU Delft Repository](https://repository.tudelft.nl/)).

## How to Implement the Model for Your Application
The main input files are located in the `grid_topology` folder. The files describing the electricity system can be adjusted for more realistic grid representations:
- `node_info.csv` --> coordinates of nodes in the system + attributed bidding zone
- `line_info.csv` --> AC lines connected between nodes with max thermal capacities
- `DC_info.csv` --> DC lines connected between nodes with max thermal capacities + From/To Direction to indicate cross-zonal interconnectors
- `incidence.csv` --> incidence matrix of all AC lines
- `incidence_dv.csv` --> incidence matrix of all DC lines
- `susceptance.csv` --> susceptance matrix of all AC lines
- `load_info.csv` --> nodal (inelastic) hourly load values
- `offshore.csv`, `onshore.csv`, and `pv.csv` --> hourly capacity factors for renewable generators per zone

Go to the `main.jl` file and adjust the loops according to the changes made. It is important to adjust the `read_data` loop based on the (created) input folder for the data. Additionally, a `Results` folder should be created where the results per simulation are stored. The default setups are `Reference_Case` (Triangular OBZ setup) and `Simple_Hybrid` (Dual OBZ setup) for the representative grid.

If a new grid topology is used, make sure to check all files included in `main.jl` for the required adjustments in the model, as indicated with comments. For computational time purposes, it is advised to simulate per month of data and copy and paste the results from the folder where results are written to a desired folder to store the results. The results can then be aggregated using `aggregate_months.py`.

A Gurobi Optimization License is required for the simulation (free for academics: [Gurobi Academic Program](https://www.gurobi.com/academia/academic-program-and-licenses/)).

After running the `main.jl` file, monthly results can be aggregated. The risk metrics can then be simulated (`analyse_risk_metrics.py`) to determine the frequency and severity of the price and volume risk. Additionally, `FB_domain.jl` provides a script for visualization of the FB domain per hour, and `create_graphs.py` provides a script for price duration curves and volume risk duration curves.

### References
M. Kenis, E. Delarue, K. Bruninx, and F. Dominguez, "Off-shore Bidding Zones Under Flow-Based Market Coupling," 2023 IEEE Belgrade PowerTech, Belgrade, Serbia, 2023, pp. 1-6, doi: 10.1109/PowerTech55446.2023.10202755.

### Contact
- LinkedIn: [Daan Verkooijen](https://www.linkedin.com/in/daan-verkooijen-9b9679189/)
- Email: dverkooijen.zakelijk@hotmail.com
