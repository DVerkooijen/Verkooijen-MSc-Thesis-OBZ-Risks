import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
import seaborn as sns

# Parameters
OBZ_zone = 4
topology = 'Reference_Case'                                             # Reference_Case, Simple_Hybrid (or new topology)
base_folder_path = 'Results\\___Case_Group1\_Transition_Mix\\total'     #Change this to select the correct folder path

if topology == 'Reference_Case':
    OWF_location = [98, 99, 103]    #Select the number in the list of nodes, so not the actual node! (python is 0 indexed, but I incroporated the 1-indexed to correspond with Julia)
                                             # For Reference_Case OWF 119, 120, 124 --> on places 98,99,103 in the nodes list 
elif topology == 'Simple_Hybrid':
    OWF_location = [98, 99]         

     # Reference_Case, _Simple_Hybrid

# Convert OWF_location to list if it is not one
if isinstance(OWF_location, int):
        OWF_location = [OWF_location]

price_df = pd.read_csv(f'{base_folder_path}/Price.csv')
res_production_df = pd.read_csv(f'{base_folder_path}/RES_production.csv')
curtailment_df = pd.read_csv(f'{base_folder_path}/Curtailment.csv')
cost_total_df = pd.read_csv(f'{base_folder_path}/redispatch/Cost_total.csv')
congestion_rent_df = pd.read_csv(f'{base_folder_path}/Congestion_Rents.csv')
net_production_OWF_df = pd.read_csv(f'{base_folder_path}/Net_Production_OWF.csv')
maxpOBZ_df = pd.read_csv(f'{base_folder_path}/FB_domain/MaxOBZ_NPs.csv')


def generate_and_save_graphs(base_folder_path, OBZ_zone):
    # Read the mean prices from the CSV file
    intra_z_indicator_pr1_df = pd.read_csv(f'{base_folder_path}/processed_results/pr1_intra_z_indicator_pr1.csv')
    mean_prices = intra_z_indicator_pr1_df['Average Price']

    # Create a color map
    colors = ['red' if i == OBZ_zone - 1 else 'black' if i == 0 else 'darkblue' if i == 1 else 'orange' for i in range(price_df.shape[0])]
    labels = ['OBZ' if i == OBZ_zone - 1 else f'Zone {i+1}' for i in range(price_df.shape[0])]

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Plot each zone's price
    for i in range(price_df.shape[0]):
        ax.plot(price_df.columns.str[1:].astype(int), price_df.iloc[i], color=colors[i], label=f'{labels[i]} (Avg: {mean_prices[i]:.2f})', linewidth=0.5, alpha=0.75, zorder=i)

        # Plot the mean price as a dotted line
        ax.axhline(mean_prices[i], color=colors[i], linestyle='dotted')

    # Set the labels
    ax.set_xlabel('Weeks')
    ax.set_ylabel('Price (€/MWh)')
    ax.set_title('Price Levels of All Zones')

    # Add a legend
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fontsize='x-small', ncol=len(labels))
    plt.tight_layout()

    # Set y-axis to start at 0
    ax.set_ylim(bottom=0, top=80)

    def hours_to_weeks(hours):
        return hours // 168 + 1

    # Set x-axis to start at 0 and set xticks to be integers
    ax.set_xlim(left=0)
    ax.set_xticks(range(0, len(price_df.columns) + 1, 168))  # Set xticks every 168 hours
    ax.set_xticklabels([hours_to_weeks(x) for x in range(0, len(price_df.columns) + 1, 168)])

    # Save the plot
    plt.savefig(f'{base_folder_path}/processed_results/price_levels.png')
    plt.close(fig)

    # For the load duration curve
    # Sort the prices in descending order
    sorted_prices = price_df.apply(lambda x: sorted(x, reverse=True), axis=1)

    # Create a new DataFrame from the sorted_prices Series
    sorted_prices_df = pd.DataFrame(sorted_prices.tolist())

    # Create a new figure and axis
    fig, ax = plt.subplots()

    # Plot each zone's sorted price
    for i in range(sorted_prices_df.shape[0]):
        ax.plot(sorted_prices_df.columns + 1, sorted_prices_df.iloc[i], color=colors[i], label=f'{labels[i]} (Avg: {mean_prices[i]:.2f})', linewidth=1.5, zorder=i)
        # Plot the mean price as a dotted line
        ax.axhline(mean_prices[i], color=colors[i], linestyle='dotted')

    # Set the labels
    ax.set_xlabel('Weeks')
    ax.set_ylabel('Price (€/MWh)')
    ax.set_title('Price Duration Curve')

    # Add a legend
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fontsize='x-small', ncol=len(labels))
    plt.tight_layout()

    # Set y-axis to start at 0
    ax.set_ylim(bottom=0, top=100)

    # Set x-axis to start at 0 and set xticks to be integers
    ax.set_xlim(left=0)
    ax.set_xticks(range(0, len(sorted_prices_df.columns) + 1, 168))  # Set xticks every 168 hours
    ax.set_xticklabels([hours_to_weeks(x) for x in range(0, len(sorted_prices_df.columns) + 1, 168)])

    # Save the plot
    plt.savefig(f'{base_folder_path}/processed_results/price_duration_curve.png')
    plt.close(fig)


def plot_volume_risk_duration_curve(base_folder_path):
    # Load the volume risk severity data
    volume_risk_severity_df = pd.read_csv(f'{base_folder_path}/processed_results/Volume_risk_Severity.csv')

    # Sort the data in descending order
    sorted_calc_risk = volume_risk_severity_df['Capacity Calculation [MWh]'].sort_values(ascending=False).reset_index(drop=True)
    sorted_alloc_risk = volume_risk_severity_df['Capacity Allocation [MWh]'].sort_values(ascending=False).reset_index(drop=True)
    sorted_total_risk = volume_risk_severity_df['Total volume risk [MWh]'].sort_values(ascending=False).reset_index(drop=True)

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Plot the data
    ax.plot(sorted_calc_risk, label='Capacity Calculation Risk [MWh]', color='red', linestyle='-', linewidth=1.5)
    ax.plot(sorted_alloc_risk, label='Capacity Allocation Risk [MWh]', color='blue', linestyle='-', linewidth=1.5)
    ax.plot(sorted_total_risk, label='Total Volume Risk [MWh]', color='green', linestyle='-', linewidth=1.5)

    # Set the labels
    ax.set_xlabel('Hours')
    ax.set_ylabel('Volume (MWh)')
    ax.set_title('Duration Curve for Volume Risks')

    # Add a legend
    ax.legend(loc='upper right')

    # Set y-axis to start at 0
    ax.set_ylim(bottom=0, top=5000)

    # Set x-axis to start at 0 and configure as needed for your data range
    ax.set_xlim(left=0, right=2000)

    #hide spines
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    # Optional: Set x-axis ticks to represent periods (weeks, months) if relevant
    # ax.set_xticks(range(0, len(sorted_total_risk), period_length))
    # ax.set_xticklabels(['Period 1', 'Period 2', 'Period 3', ...])

    plt.grid(True)
    plt.tight_layout()

    # Save the plot
    plt.savefig(f'{base_folder_path}/processed_results/volume_risk_duration_curve.png')

    plt.close(fig)


plot_volume_risk_duration_curve(base_folder_path)

generate_and_save_graphs(base_folder_path, OBZ_zone)
