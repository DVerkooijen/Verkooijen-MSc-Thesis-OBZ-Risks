import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os

# Parameters
OBZ_zone = 4
topology = 'Simple_Hybrid'                                       # Reference_Case, Simple_Hybrid (or new topology)
base_folder_path = 'Results\\___Case_Group3\\Asym. Z1\\total'    #Change this to select the correct folder path
hydrogen = "yes"   
WTP = 41.70                                                      #Change this to the correct WTP value (in case of demand agent) 

if topology == 'Reference_Case':
    OWF_location = [98, 99, 103]    #Select the number in the list of nodes, so not the actual node! (python is 0 indexed, but I incroporated the 1-indexed to correspond with Julia)
                                             # For Reference_Case OWF 119, 120, 124 --> on places 98,99,103 in the nodes list 
elif topology == 'Simple_Hybrid':
    OWF_location = [98, 99]         

# Convert OWF_location to list if it is not one
if isinstance(OWF_location, int):
        OWF_location = [OWF_location]

# Read all the main data files
price_df = pd.read_csv(f'{base_folder_path}/Price.csv')
res_production_df = pd.read_csv(f'{base_folder_path}/RES_production.csv')
curtailment_df = pd.read_csv(f'{base_folder_path}/Curtailment.csv')
cost_total_df = pd.read_csv(f'{base_folder_path}/redispatch/Cost_total.csv')
congestion_rent_df = pd.read_csv(f'{base_folder_path}/Congestion_Rents.csv')
net_production_OWF_df = pd.read_csv(f'{base_folder_path}/Net_Production_OWF.csv')
maxpOBZ_df = pd.read_csv(f'{base_folder_path}/FB_domain/MaxOBZ_NPs.csv')


def compute_general_results(OWF_location, base_folder_path):

    # Revenue calculation considering curtailment
    price_OBZ = price_df.iloc[OBZ_zone - 1].values
    total_revenue_per_hour = pd.Series([0] * len(price_OBZ))
   
    for node in OWF_location:
        net_production_OWF = net_production_OWF_df.iloc[node-1].values
        total_revenue_per_hour += price_OBZ * net_production_OWF

    revenue_df = pd.DataFrame({
        'Hour': range(1, len(total_revenue_per_hour) + 1),
        'Revenue': total_revenue_per_hour
    })
    total_revenue = total_revenue_per_hour.sum()

    # Additional calculations
    total_curtailment = curtailment_df.iloc[[loc-1 for loc in OWF_location]].sum().sum()
    total_day_ahead_cost = cost_total_df['cost_generation'].sum()
    total_redispatch_cost = abs(cost_total_df['cost_redispatch'].sum())
    # total_cost =cost_total_df['total_cost'].sum()
    total_congestion_rent = congestion_rent_df['Rent'].abs().sum()

    # Calculate wind metrics at zero and positive prices
    zero_price_hours = np.isclose(price_OBZ, 0, atol=1e-4)
    positive_price_hours = price_OBZ > 0
    wind_at_zero_price = sum(sum(net_production_OWF_df.iloc[loc-1].values * zero_price_hours) for loc in OWF_location)
    wind_at_positive_price = sum(sum(net_production_OWF_df.iloc[loc-1].values * positive_price_hours) for loc in OWF_location)
    total_production = sum(net_production_OWF_df.iloc[:, loc-1].sum() for loc in OWF_location)
    available_wind = sum(res_production_df.iloc[:, loc-1].sum() for loc in OWF_location)

    # Consolidate results
    results = {
        'Metrics': ['Total Revenue [M€]','Available Wind [GWh]', 'Total Production [GWh]', 'Total Curtailment [GWh]', 'Wind at Zero-Price [GWh]', 'Wind at Positive Price [GWh]',
                    'Total Congestion Rent [M€]', 'Total Day-Ahead Cost [M€]', 'Total Redispatch Cost [M€]'],
        'Values': [total_revenue / 1_000_000,  # Convert to M€
                available_wind/1_000,  # Convert to GWh   
                total_production/1_000,  # Convert to GWh
                total_curtailment / 1_000,  # Convert to GWh
                wind_at_zero_price / 1_000,  # Convert to GWh
                wind_at_positive_price / 1_000,  # Convert to GWh
                total_congestion_rent / 1_000_000,  # Convert to M€
                total_day_ahead_cost / 1_000_000,  # Convert to M€
                total_redispatch_cost / 1_000_000]  # Convert to M€
    }
    results_df = pd.DataFrame(results)

    # # Print results to the console
    # print(revenue_df)
    # print("Total Revenue: €", total_revenue)
    # print("Total Day-Ahead Cost: €", total_day_ahead_cost)
    # print("Total Redispatch Cost: €", total_redispatch_cost)
    # print("Total Congestion Rent: €", total_congestion_rent)
    # print("Total Wind at Zero-Price: MWh", wind_at_zero_price)
    # print("Total Wind at Price > 0 for nodes", OWF_location, ": MWh", wind_at_positive_price)
    # print(results_df)

    results_df.to_csv(f'{base_folder_path}/processed_results/_General_Results.csv', index=False)


#Price Risk 1: Flat Price risk
def compute_price_risk_1(OBZ_zone, base_folder_path):
    mean_price = price_df.mean(axis=1)
    std_dev_price = price_df.std(axis=1)
    coefficient_of_variation = std_dev_price / mean_price

    intra_z_indicator_pr1_df = pd.DataFrame({
        'Zone': range(1, len(mean_price) + 1),
        'Average Price': mean_price.values,
        'Std Deviation': std_dev_price.values,
        'Coefficient of Variation': coefficient_of_variation.values
    })
    intra_z_indicator_pr1_df['Zone'] = intra_z_indicator_pr1_df['Zone'].apply(lambda x: 'OBZ' if x == OBZ_zone else x)
    intra_z_indicator_pr1_df.to_csv(f'{base_folder_path}/processed_results/pr1_intra_z_indicator_pr1.csv', index=False)

    inter_zone_mean = mean_price.mean()
    inter_zone_std_dev = mean_price.std()
    inter_zone_coefficient_of_variation = inter_zone_std_dev / inter_zone_mean

    # Create DataFrame for inter-zone indicators
    inter_z_indicator_pr1_df = pd.DataFrame({
        'Metric': ['Inter-Zone Mean Price', 'Inter-Zone Std Deviation', 'Inter-Zone Coefficient of Variation'],
        'Value': [inter_zone_mean, inter_zone_std_dev, inter_zone_coefficient_of_variation]
    })

    # Save to CSV
    inter_z_indicator_pr1_df.to_csv(f'{base_folder_path}/processed_results/pr1_inter_z_indicator.csv', index=False)

#price risk 2
def compute_price_risk_2(OBZ_zone, base_folder_path):
    # Initialize counters
    positive_hours = 0
    negative_hours = 0
    positive_mwh = 0
    negative_mwh = 0
    positive_effect = 0
    negative_effect = 0

    # Iterate over each hour
    for i in range(len(price_df.columns)):
        price_OBZ = np.round(price_df.iloc[OBZ_zone - 1, i], 4)
        prices_other_zones = np.round(price_df.drop(OBZ_zone - 1).iloc[:, i], 4)

        #skip hydrogen WTP hours
        if price_OBZ == WTP:
            continue

        # Check if price OBZ != price other zones & also not 0
        if not prices_other_zones.eq(price_OBZ).any() and abs(round(price_OBZ, 4)) != 0:
            # Calculate net production for this hour
            net_production = sum(net_production_OWF_df.iloc[loc-1, i] for loc in OWF_location)

            # Check if it's a positive or negative hour
            if abs(round(price_OBZ, 4)) > prices_other_zones.min():
                positive_hours += 1
                positive_mwh += net_production
                positive_effect += net_production * abs(round(price_OBZ, 4))
            elif abs(round(price_OBZ, 4)) < prices_other_zones.min():
                negative_hours += 1
                negative_mwh += net_production
                negative_effect += net_production * abs(round(price_OBZ, 4))

    # Calculate net effect
    net_effect = positive_effect - negative_effect

    # Consolidate results
    results = {
        'Metrics': ['Number of Positive Hours', 'Number of Negative Hours', 'Positive Hours [MWh]', 'Negative Hours [MWh]',
                    'Positive Effect [€]', 'Negative Effect [€]', 'Net Effect [€]'],
        'Values': [positive_hours, negative_hours, positive_mwh, negative_mwh, positive_effect, negative_effect, net_effect]
    }
    results_df = pd.DataFrame(results)

    # Save to CSV
    results_df.to_csv(f'{base_folder_path}/processed_results/pr2_results.csv', index=False)

if topology == 'Simple_Hybrid':
    def compute_price_risk_3(OBZ_zone, OWF_location, base_folder_path):
        # Initialize counters
        zero_price_hours = 0
        zero_price_mwh = 0
        zero_price_hour_RES = 0
        zero_price_mwh_RES = 0
        zero_price_curt_mwh = 0
        zero_price_curt_mwh_RES = 0
        price_convergence_z1_hours = 0 
        price_convergence_z1_mwh = 0  
        price_convergence_z2_hours = 0 
        price_convergence_z2_mwh = 0
        price_convergence_z3_hours = 0
        price_convergence_z3_mwh = 0
        other_price_formation_hours = 0
        other_price_formation_mwh = 0

        hourly_data = []
        # Iterate over each hour
        for i in range(len(price_df.columns)):
            price_OBZ = price_df.iloc[OBZ_zone - 1, i]
            price_z1 = price_df.iloc[0, i]
            price_z2 = price_df.iloc[1, i]
            price_z3 = price_df.iloc[2, i]
            prices_other_zones = price_df.drop(OBZ_zone - 1).iloc[:, i]

            # Check zero-price conditions due to RES
            if abs(round(price_OBZ, 4)) == 0 and prices_other_zones.apply(lambda x: round(x, 4)).eq(0).any():
                zero_price_hour_RES += 1
                zero_price_mwh_RES += sum(net_production_OWF_df.iloc[loc-1, i] for loc in OWF_location)
                zero_price_curt_mwh_RES += sum(curtailment_df.iloc[loc-1, i] for loc in OWF_location)
                hourly_data.append({'Hour': i+1, 'Category': 'Zero-Price RES', 'Curtailed MWh': sum(curtailment_df.iloc[loc-1, i] for loc in OWF_location)})

            # Check zero-price conditions not due to RES
            elif abs(round(price_OBZ, 4)) == 0 and not prices_other_zones.apply(lambda x: round(x, 4)).eq(0).any():
                zero_price_hours += 1
                zero_price_mwh += sum(net_production_OWF_df.iloc[loc-1, i] for loc in OWF_location)
                zero_price_curt_mwh += sum(curtailment_df.iloc[loc-1, i] for loc in OWF_location)
                hourly_data.append({'Hour': i+1, 'Category': 'Zero-Price', 'Curtailed MWh': sum(curtailment_df.iloc[loc-1, i] for loc in OWF_location)})

            # Check price convergence and other price formations
            else:
                if round(price_OBZ, 2) == round(price_z1, 2):
                        price_convergence_z1_hours += 1
                        price_convergence_z1_mwh += sum(net_production_OWF_df.iloc[loc-1, i] for loc in OWF_location)
                        hourly_data.append({'Hour': i+1, 'Category': 'Price Convergence Zone 1', 'Curtailed MWh': sum(curtailment_df.iloc[loc-1, i] for loc in OWF_location)})
                elif round(price_OBZ, 2) == round(price_z2, 2):
                        price_convergence_z2_hours += 1
                        price_convergence_z2_mwh += sum(net_production_OWF_df.iloc[loc-1, i] for loc in OWF_location)
                        hourly_data.append({'Hour': i+1, 'Category': 'Price Convergence Zone 2', 'Curtailed MWh': sum(curtailment_df.iloc[loc-1, i] for loc in OWF_location)})
                elif round(price_OBZ, 2) == round(price_z3, 2):
                        price_convergence_z3_hours += 1
                        price_convergence_z3_mwh += sum(net_production_OWF_df.iloc[loc-1, i] for loc in OWF_location)
                        hourly_data.append({'Hour': i+1, 'Category': 'Price Convergence Zone 3', 'Curtailed MWh': sum(curtailment_df.iloc[loc-1, i] for loc in OWF_location)})
                else:
                        other_price_formation_hours += 1
                        other_price_formation_mwh += sum(net_production_OWF_df.iloc[loc-1, i] for loc in OWF_location)
                        hourly_data.append({'Hour': i+1, 'Category': 'Other Price Formation', 'Curtailed MWh': sum(curtailment_df.iloc[loc-1, i] for loc in OWF_location)})

        # Consolidate results
        results = {
                'Metrics': [
                    '#Price Collapse Hours', 'Price Collapse Hours [MWh]', 'Curtailment [MWh]',
                    '#Zero-Priced Hours due to RES', 'Zero-Priced Hours due to RES [MWh]', 'Curtailment due to RES [MWh]',
                    '#Price Convergence Zone 1 Hours', 'Price Convergence Zone 1 MWh',
                    '#Price Convergence Zone 2 Hours', 'Price Convergence Zone 2 MWh',
                    '#Price Convergence Zone 3 Hours', 'Price Convergence Zone 3 MWh',
                    '#Other Price Formation Hours', 'Other Price Formation MWh'
                ],
                'Values': [
                    zero_price_hours, zero_price_mwh, zero_price_curt_mwh,
                    zero_price_hour_RES, zero_price_mwh_RES, zero_price_curt_mwh_RES,
                    price_convergence_z1_hours, price_convergence_z1_mwh,
                    price_convergence_z2_hours, price_convergence_z2_mwh,
                    price_convergence_z3_hours, price_convergence_z3_mwh,
                    other_price_formation_hours, other_price_formation_mwh
                ]
            }
        results_df = pd.DataFrame(results)

        os.makedirs(f'{base_folder_path}/processed_results', exist_ok=True)
        results_df.to_csv(f'{base_folder_path}/processed_results/pr3_results.csv', index=False)

        hourly_data_df = pd.DataFrame(hourly_data)
        hourly_data_df.to_csv(f'{base_folder_path}/processed_results/pr3_results_hourly.csv', index=False)

if topology == 'Reference_Case':
    def compute_price_risk_3(OBZ_zone, OWF_location, base_folder_path):
        # Initialize counters
        zero_price_hours = 0
        zero_price_mwh = 0
        zero_price_hour_RES = 0
        zero_price_mwh_RES = 0
        zero_price_curt_mwh = 0
        zero_price_curt_mwh_RES = 0
        price_convergence_low_hours = 0 
        price_convergence_low_mwh = 0  
        price_convergence_high_hours = 0
        price_convergence_high_mwh = 0
        price_convergence_middle_hours = 0
        price_convergence_middle_mwh = 0
        other_price_formation_hours = 0
        other_price_formation_mwh = 0

        hourly_data = []
        # Iterate over each hour
        for i in range(len(price_df.columns)):
            price_OBZ = price_df.iloc[OBZ_zone - 1, i]
            prices_other_zones = price_df.drop(OBZ_zone - 1).iloc[:, i]
        
            # Check zero-price conditions due to RES
            if abs(round(price_OBZ, 4)) == 0 and prices_other_zones.apply(lambda x: round(x, 4)).eq(0).any():
                zero_price_hour_RES += 1
                zero_price_mwh_RES += sum(net_production_OWF_df.iloc[loc-1, i] for loc in OWF_location)
                zero_price_curt_mwh_RES += sum(curtailment_df.iloc[loc-1, i] for loc in OWF_location)
                hourly_data.append({'Hour': i+1, 'Category': 'Zero-Price RES', 'Curtailed MWh': sum(curtailment_df.iloc[loc-1, i] for loc in OWF_location)})

            # Check zero-price conditions not due to RES
            elif abs(round(price_OBZ, 4)) == 0 and not prices_other_zones.apply(lambda x: round(x, 4)).eq(0).any():
                zero_price_hours += 1
                zero_price_mwh += sum(net_production_OWF_df.iloc[loc-1, i] for loc in OWF_location)
                zero_price_curt_mwh += sum(curtailment_df.iloc[loc-1, i] for loc in OWF_location)
                hourly_data.append({'Hour': i+1, 'Category': 'Zero-Price', 'Curtailed MWh': sum(curtailment_df.iloc[loc-1, i] for loc in OWF_location)})

            # Check price convergence and other price formations
            else:
                min_price = prices_other_zones.min()
                max_price = prices_other_zones.max()
                other_prices = [x for x in prices_other_zones if x != min_price and x != max_price]

                if round(price_OBZ, 2) == round(min_price, 2):
                    price_convergence_low_hours += 1
                    price_convergence_low_mwh += sum(net_production_OWF_df.iloc[loc-1, i] for loc in OWF_location)
                    hourly_data.append({'Hour': i+1, 'Category': 'Price Convergence Low', 'Curtailed MWh': sum(curtailment_df.iloc[loc-1, i] for loc in OWF_location)})
                elif round(price_OBZ, 2) == round(max_price, 2):
                    price_convergence_high_hours += 1
                    price_convergence_high_mwh += sum(net_production_OWF_df.iloc[loc-1, i] for loc in OWF_location)
                    hourly_data.append({'Hour': i+1, 'Category': 'Price Convergence High', 'Curtailed MWh': sum(curtailment_df.iloc[loc-1, i] for loc in OWF_location)})
                elif any(round(price_OBZ, 2) == round(x, 2) for x in other_prices):
                    price_convergence_middle_hours += 1
                    price_convergence_middle_mwh += sum(net_production_OWF_df.iloc[loc-1, i] for loc in OWF_location)
                    hourly_data.append({'Hour': i+1, 'Category': 'Price Convergence Middle', 'Curtailed MWh': sum(curtailment_df.iloc[loc-1, i] for loc in OWF_location)})
                else:
                    other_price_formation_hours += 1
                    other_price_formation_mwh += sum(net_production_OWF_df.iloc[loc-1, i] for loc in OWF_location)
                    hourly_data.append({'Hour': i+1, 'Category': 'Other Price Formation', 'Curtailed MWh': sum(curtailment_df.iloc[loc-1, i] for loc in OWF_location)})

        # Consolidate results
        results = {
            'Metrics': [
                '#Price Collapse Hours', 'Price Collapse Hours [MWh]', 'Curtailment [MWh]',
                '#Zero-Priced Hours due to RES', 'Zero-Priced Hours due to RES [MWh]', 'Curtailment due to RES [MWh]',
                '#Price Convergence Low Hours', 'Price Convergence Low MWh',
                '#Price Convergence High Hours', 'Price Convergence High MWh',
                '#Price Convergence Middle Hours', 'Price Convergence Middle MWh',
                '#Other Price Formation Hours', 'Other Price Formation MWh'
            ],
            'Values': [
                zero_price_hours, zero_price_mwh, zero_price_curt_mwh,
                zero_price_hour_RES, zero_price_mwh_RES, zero_price_curt_mwh_RES,
                price_convergence_low_hours, price_convergence_low_mwh,
                price_convergence_high_hours, price_convergence_high_mwh,
                price_convergence_middle_hours, price_convergence_middle_mwh,
                other_price_formation_hours, other_price_formation_mwh
            ]
        }
        results_df = pd.DataFrame(results)

        os.makedirs(f'{base_folder_path}/processed_results', exist_ok=True)
        results_df.to_csv(f'{base_folder_path}/processed_results/pr3_results.csv', index=False)

        hourly_data_df = pd.DataFrame(hourly_data)
        hourly_data_df.to_csv(f'{base_folder_path}/processed_results/pr3_results_hourly.csv', index=False)



def compute_volume_risks(base_folder_path, OBZ_zone, OWF_location):
    # Initialize data frames to store results
    volume_risk_severity = []
    volume_risk_frequency = {
        'total': 0,
        'Capacity Calculation': 0,
        'Capacity Allocation': 0,
        'both': 0
    }

    # Loop through each hour
    for hour in price_df.columns:
        maxpOBZ = maxpOBZ_df[hour].iloc[OBZ_zone - 1]
        potential_RES = sum(res_production_df[hour].iloc[owf - 1] for owf in OWF_location)
        actual_RES = sum(net_production_OWF_df[hour].iloc[owf - 1] for owf in OWF_location)
        curtailment = sum(curtailment_df[hour].iloc[owf - 1] for owf in OWF_location)
        price = price_df[hour].iloc[OBZ_zone - 1]

        # Get the prices for each zone
        price_zone1 = price_df[hour].iloc[0]
        price_zone2 = price_df[hour].iloc[1]
        price_zone3 = price_df[hour].iloc[2]

        if any(curtailment_df[hour].iloc[owf - 1] != 0 for owf in OWF_location):
            # Calculate the amount of MWs curtailed for each wind farm
            curtailment_per_owf = [curtailment_df[hour].iloc[owf - 1] for owf in OWF_location]

            if topology == 'Reference_Case':	
                if maxpOBZ < potential_RES and maxpOBZ >= actual_RES:
                    V3 = potential_RES - maxpOBZ
                    V4 = curtailment - V3

                    if V3 > 0 and V4 == 0:
                        volume_risk_frequency['Capacity Calculation'] += 1
                        volume_risk_frequency['total'] += 1
                        volume_risk_severity.append({
                            'Hour': hour,
                            'Capacity Calculation [MWh]': V3,
                            'Capacity Allocation [MWh]': 0,
                            'Total volume risk [MWh]': V3,
                            'λ OBZ [€/MWh]': price,
                            'λ zone 1 [€/MWh]': price_zone1,
                            'λ zone 2 [€/MWh]': price_zone2,
                            'λ zone 3 [€/MWh]': price_zone3,
                            'Curtailment OWF 119 (z1) [MWh]': curtailment_per_owf[0],
                            'Curtailment OWF 124 (z2) [MWh]': curtailment_per_owf[2],
                            'Curtailment OWF 120 (z3) [MWh]': curtailment_per_owf[1]
                        })
                    elif V3 > 0 and V4 > 0:
                        volume_risk_frequency['both'] += 1
                        volume_risk_frequency['Capacity Calculation'] += 1
                        volume_risk_frequency['Capacity Allocation'] += 1
                        volume_risk_frequency['total'] += 1
                        volume_risk_severity.append({
                            'Hour': hour,
                            'Capacity Calculation [MWh]': V3,
                            'Capacity Allocation [MWh]': V4,
                            'Total volume risk [MWh]': V3 + V4,
                            'λ OBZ [€/MWh]': price,
                            'λ zone 1 [€/MWh]': price_zone1,
                            'λ zone 2 [€/MWh]': price_zone2,
                            'λ zone 3 [€/MWh]': price_zone3,
                            'Curtailment OWF 119 (z1) [MWh]': curtailment_per_owf[0],
                            'Curtailment OWF 124 (z2) [MWh]': curtailment_per_owf[2],
                            'Curtailment OWF 120 (z3) [MWh]': curtailment_per_owf[1]
                        })
                elif maxpOBZ >= potential_RES and maxpOBZ >= actual_RES:
                    V4 = curtailment
                    V3 = 0
                    volume_risk_frequency['Capacity Allocation'] += 1
                    volume_risk_frequency['total'] += 1
                    volume_risk_severity.append({
                        'Hour': hour,
                        'Capacity Calculation [MWh]': 0,
                        'Capacity Allocation [MWh]': V4,
                        'Total volume risk [MWh]': V4,
                        'λ OBZ [€/MWh]': price,
                        'λ zone 1 [€/MWh]': price_zone1,
                        'λ zone 2 [€/MWh]': price_zone2,
                        'λ zone 3 [€/MWh]': price_zone3,
                        'Curtailment OWF 119 (z1) [MWh]': curtailment_per_owf[0],
                        'Curtailment OWF 124 (z2) [MWh]': curtailment_per_owf[2],
                        'Curtailment OWF 120 (z3) [MWh]': curtailment_per_owf[1]
                    })

            elif topology == 'Simple_Hybrid':
                if maxpOBZ < potential_RES and maxpOBZ >= actual_RES:
                    V3 = potential_RES - maxpOBZ
                    V4 = curtailment - V3

                    if V3 > 0 and V4 == 0:
                        volume_risk_frequency['Capacity Calculation'] += 1
                        volume_risk_frequency['total'] += 1
                        volume_risk_severity.append({
                            'Hour': hour,
                            'Capacity Calculation [MWh]': V3,
                            'Capacity Allocation [MWh]': 0,
                            'Total volume risk [MWh]': V3,
                            'λ OBZ [€/MWh]': price,
                            'λ zone 1 [€/MWh]': price_zone1,
                            'λ zone 2 [€/MWh]': price_zone2,
                            'λ zone 3 [€/MWh]': price_zone3,
                            'Curtailment OWF 119 (z1) [MWh]': curtailment_per_owf[0],
                            'Curtailment OWF 120 (z3) [MWh]': curtailment_per_owf[1]
                        })
                    elif V3 > 0 and V4 > 0:
                        volume_risk_frequency['both'] += 1
                        volume_risk_frequency['Capacity Calculation'] += 1
                        volume_risk_frequency['Capacity Allocation'] += 1
                        volume_risk_frequency['total'] += 1
                        volume_risk_severity.append({
                            'Hour': hour,
                            'Capacity Calculation [MWh]': V3,
                            'Capacity Allocation [MWh]': V4,
                            'Total volume risk [MWh]': V3 + V4,
                            'λ OBZ [€/MWh]': price,
                            'λ zone 1 [€/MWh]': price_zone1,
                            'λ zone 2 [€/MWh]': price_zone2,
                            'λ zone 3 [€/MWh]': price_zone3,
                            'Curtailment OWF 119 (z1) [MWh]': curtailment_per_owf[0],
                            'Curtailment OWF 120 (z3) [MWh]': curtailment_per_owf[1]
                        })
                elif maxpOBZ >= potential_RES and maxpOBZ >= actual_RES:
                    V4 = curtailment
                    V3 = 0
                    volume_risk_frequency['Capacity Allocation'] += 1
                    volume_risk_frequency['total'] += 1
                    volume_risk_severity.append({
                        'Hour': hour,
                        'Capacity Calculation [MWh]': 0,
                        'Capacity Allocation [MWh]': V4,
                        'Total volume risk [MWh]': V4,
                        'λ OBZ [€/MWh]': price,
                        'λ zone 1 [€/MWh]': price_zone1,
                        'λ zone 2 [€/MWh]': price_zone2,
                        'λ zone 3 [€/MWh]': price_zone3,
                        'Curtailment OWF 119 (z1) [MWh]': curtailment_per_owf[0],
                        'Curtailment OWF 120 (z3) [MWh]': curtailment_per_owf[1]
                    })

    # Convert results to DataFrame
    volume_risk_severity_df = pd.DataFrame(volume_risk_severity)
    volume_risk_frequency_df = pd.DataFrame([volume_risk_frequency])

    # Save results to CSV
    volume_risk_severity_df.to_csv(os.path.join(base_folder_path, 'processed_results/Volume_risk_Severity.csv'), index=False)
    volume_risk_frequency_df.to_csv(os.path.join(base_folder_path, 'processed_results/Volume_risk_Frequency.csv'), index=False)




# # Call the function
compute_general_results(OWF_location, base_folder_path)

compute_price_risk_1(OBZ_zone, base_folder_path)
compute_price_risk_2(OBZ_zone, base_folder_path)
compute_price_risk_3(OBZ_zone, OWF_location, base_folder_path)

# Execute the function
compute_volume_risks(base_folder_path, OBZ_zone, OWF_location)




