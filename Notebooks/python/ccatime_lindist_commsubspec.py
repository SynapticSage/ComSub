import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.utils import resample
from sklearn.preprocessing import MinMaxScaler
from tqdm import tqdm

folder = '/Volumes/MATLAB-Drive/Shared/figures/tables/'
name   = 'ZT2powerccatime'
df_csv = pd.read_csv(os.path.join(folder, f'{name}.csv'))

# Number of bins and bootstrap samples
n_bins = 100
n_bootstrap_samples = 1000

# Update the columns to bootstrap
columns_to_bootstrap = ["U1", "U2", "U3", "V1", "V2", "V3", "Cavgtheta", 
                        "Cavgdelta", "Cavgripple", "S1theta", "S1delta", 
                        "S2theta", "wpli_avgtheta", "S1ripple", "S2ripple",
                        "S2delta", "wpli_avgripple", "wpli_avgdelta"]

# Reinitialize the DataFrame to hold the results
bootstrap_means_combined = []

# Rerun the bootstrap for each bin and each column, for each trajectory bound
for trajbound in tqdm([0, 1], desc="trajbound", total=2):
    # Filter the data based on trajbound
    data_trajbound = df_csv[df_csv["trajbound"] == trajbound]
    # Bin the lindist data
    data_trajbound["lindist_bin"] = pd.cut(data_trajbound["lindist"], bins=n_bins)
    
    for bin_label in tqdm(data_trajbound["lindist_bin"].unique().categories, desc="lindist_bin", total=n_bins):
        for column in columns_to_bootstrap:
            # Get the data for this bin and column
            data = data_trajbound.loc[data_trajbound["lindist_bin"] == bin_label, column]
            # Generate bootstrap samples and compute their means
            bootstrap_means = [resample(data).mean() for _ in range(n_bootstrap_samples)]
            # Add the results to the DataFrame
            for bootstrap_mean in bootstrap_means:
                bootstrap_means_combined.append({
                    "lindist_bin": bin_label,
                    "column": column,
                    "bootstrap_mean": bootstrap_mean,
                    "trajbound": trajbound
                })

# Convert the list of dictionaries to a DataFrame
bootstrap_means_combined = pd.DataFrame(bootstrap_means_combined)
bootstrap_means_combined.to_csv(os.path.join(folder, f'{name}_bootstrap.csv'), index=False)
# ----------------------------------------------------

# Normalize the bootstrap_mean values
scaler = MinMaxScaler()
for column in tqdm(columns_to_bootstrap, desc="feature engineering", total=len(columns_to_bootstrap)):
    # Get the data for this column
    data = bootstrap_means_combined[bootstrap_means_combined["column"] == column]["bootstrap_mean"].values.reshape(-1, 1)
    # Scale the data
    scaled_data = scaler.fit_transform(data)
    # Update the DataFrame
    bootstrap_means_combined.loc[bootstrap_means_combined["column"] == column, "bootstrap_mean"] = scaled_data

bootstrap_means_combined.head()
bootstrap_means_combined.to_csv(os.path.join(folder, f'{name}_bootstrap_normalized.csv'), index=False)
# ----------------------------------------------------


# Update the components for each row of the subplot grid
row_components = [["U1", "U2", "U3"], ["V1", "V2", "V3"], ["Cavgtheta", "S1theta", "S2theta", "wpli_avgtheta"],
                  ["Cavgdelta", "S1delta", "S2delta", "wpli_avgdelta"], ["Cavgripple", "S1ripple", "S2ripple", "wpli_avgripple"]]


# Add the lindist_bin_mid column back to the DataFrame
bootstrap_means_combined["lindist_bin_mid"] = bootstrap_means_combined["lindist_bin"].apply(lambda x: x.mid)

# Create a 5x2 grid of subplots
fig, axes = plt.subplots(nrows=5, ncols=2, figsize=(15, 25))

# Set the overall title
fig.suptitle('Normalized Bootstrap Means for Different Components and Trajectories', fontsize=20)

for i, components in enumerate(row_components):
    for j, trajbound in enumerate(column_trajbounds):
        # Get the data for this subplot
        data = bootstrap_means_combined[
            bootstrap_means_combined["column"].isin(components) &
            (bootstrap_means_combined["trajbound"] == trajbound)
        ]
        # Create the subplot
        sns.lineplot(x="lindist_bin_mid", y="bootstrap_mean", hue="column", data=data, ax=axes[i, j])
        # Set the title and labels
        axes[i, j].set_title(f'Trajbound = {trajbound}')
        axes[i, j].set_xlabel("Lindist Bin Midpoint")
        axes[i, j].set_ylabel("Normalized Bootstrap Mean")
        # Rotate the x-axis labels for readability
        axes[i, j].tick_params(axis='x', rotation=45)

# Improve the layout
plt.tight_layout()
# Add space for the overall title
fig.subplots_adjust(top=0.92)
plt.show()

# ----------------------------------------------------

# Set the flag for shading the confidence intervals to False
shade_confidence_intervals = False
stratify_spectral_components = True  # Change this to False if you don't want to stratify and shade the spectral components

# Create a 5x2 grid of subplots
fig, axes = plt.subplots(nrows=5, ncols=2, figsize=(15, 25))

# Set the overall title
fig.suptitle('Normalized Bootstrap Means for Different Components and Trajectories', fontsize=20)

for i, components in enumerate(row_components):
    for j, trajbound in enumerate(column_trajbounds):
        # Get the data for this subplot
        data = bootstrap_means_combined[
            bootstrap_means_combined["column"].isin(components) &
            (bootstrap_means_combined["trajbound"] == trajbound)
        ]
        # Create the subplot
        for component in components:
            component_data = data[data["column"] == component].sort_values(by="lindist_bin_mid")
            # Add the base value to the bootstrap_mean if stratification is active
            if stratify_spectral_components and component.startswith(("S1", "S2")):
                component_data["bootstrap_mean"] += component_fill_bases[component]
            # Plot the curve
            axes[i, j].plot(component_data["lindist_bin_mid"], component_data["bootstrap_mean"], color=component_colors[component], label=component)
            # Fill under the curve based on the component type and the flag
            if component in ["U1", "U2", "U3", "V1", "V2", "V3"] or (stratify_spectral_components and component.startswith(("S1", "S2"))):
                axes[i, j].fill_between(component_data["lindist_bin_mid"], component_fill_bases[component], component_data["bootstrap_mean"], color=component_colors[component], alpha=0.3)
            # Shade the confidence interval around the curve if the flag is set
            if shade_confidence_intervals:
                component_data_bootstrap_samples = component_data["bootstrap_mean"].tolist()
                confidence_interval = np.percentile(component_data_bootstrap_samples, [2.5, 97.5])
                axes[i, j].fill_between(component_data["lindist_bin_mid"], confidence_interval[0], confidence_interval[1], color=component_colors[component], alpha=0.1)
        # Set the title and labels
        axes[i, j].set_title(f'Trajbound = {trajbound}')
        axes[i, j].set_xlabel("Lindist Bin Midpoint")
        axes[i, j].set_ylabel("Normalized Bootstrap Mean")
        # Rotate the x-axis labels for readability
        axes[i, j].tick_params(axis='x', rotation=45)
        # Add a legend
        axes[i, j].legend()

# Improve the layout
plt.tight_layout()
# Add space for the overall title
fig.subplots_adjust(top=0.92)
plt.show()
