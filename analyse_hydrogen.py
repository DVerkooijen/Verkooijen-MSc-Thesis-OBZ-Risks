import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

base_folder_path = 'Results\\___Case_Group5\\_Off_1WTP43.77\\total'      #Change this to select the correct folder path


OWF_location = [97, 98, 102]  #Select the number in the list of nodes, so not the actual node! For OWF 119, 120, 124 --> on nodes 98,99,103
OBZ_zone = 3

WTP = 43.77

price_df = pd.read_csv(f'{base_folder_path}/Price.csv')
res_production_df = pd.read_csv(f'{base_folder_path}/RES_production.csv')
curtailment_df = pd.read_csv(f'{base_folder_path}/Curtailment.csv')
cost_total_df = pd.read_csv(f'{base_folder_path}/redispatch/Cost_total.csv')
congestion_rent_df = pd.read_csv(f'{base_folder_path}/Congestion_Rents.csv')
net_production_OWF_df = pd.read_csv(f'{base_folder_path}/Net_Production_OWF.csv')
maxpOBZ_df = pd.read_csv(f'{base_folder_path}/FB_domain/MaxOBZ_NPs.csv')
h2_dispatch_df = pd.read_csv(f'{base_folder_path}/Hydrogen_dispatch.csv')

h2_hours = price_df.iloc[:, OBZ_zone].eq(WTP).sum()
df_to_save = pd.DataFrame({'H2_Hours': [h2_hours]})
# df_to_save.to_csv('h2_hours.csv', index=False)

print(f'Hours with hydrogen price of {WTP} €/MWh: {h2_hours}')


FLH = []
NFLH = []
ND = []
for index, row in h2_dispatch_df.iterrows():
    FLH.append((row == 1.0).sum())
    NFLH.append(((row > 0) & (row < 1)).sum())
    ND.append((row == 0).sum())

Operation_df = pd.DataFrame({'FLH': FLH, 'NFLH': NFLH, 'ND': ND	})
Operation_df.to_csv(f'{base_folder_path}/processed_results/H2_Operation.csv', index=False)


h2_hours_equal_WTP = 0
total_exported_power_during_h2_hours = 0

# Iterate through each column in price_df (excluding the first column if it's not a price column)
for col in range(price_df.shape[1]):
    if price_df.iloc[OBZ_zone, col] == WTP:
        h2_hours_equal_WTP += 1
        # Sum the production from the specified OWFs for this hour
        exported_power_this_hour = net_production_OWF_df.iloc[OWF_location, col].sum()
        total_exported_power_during_h2_hours += exported_power_this_hour

total_exported_power_during_h2_hours /= 1000

# Prepare and save results
data = {
    'H2 Hours Equal WTP': [h2_hours_equal_WTP],
    'Total Export H2 Hours': [total_exported_power_during_h2_hours]
}
results_df = pd.DataFrame(data)
results_df.to_csv(f'{base_folder_path}/processed_results/H2_prices.csv', index=False)
