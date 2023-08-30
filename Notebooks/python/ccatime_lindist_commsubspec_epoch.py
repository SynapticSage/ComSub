import os, gc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
from sklearn.utils import resample
from sklearn.preprocessing import MinMaxScaler
from tqdm import tqdm
import itertools
import matplotlib.cm as cm

def downcast_dataframe(df):
    """
    Downcasts a DataFrame's columns to save memory.
    
    - Int64/Int32 columns are casted to Int16.
    - Float64 columns are casted to Float32.

    Args:
    - df (pd.DataFrame): The input dataframe.

    Returns:
    - pd.DataFrame: The downcasted dataframe.
    """
    
    for col in df.columns:
        if df[col].dtype == 'int64' or df[col].dtype == 'int32':
            df[col] = df[col].astype('int16')
        elif df[col].dtype == 'float64':
            df[col] = df[col].astype('float32')
    
    return df


# Set the flag for shading the confidence intervals to False
# - - - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - - 
# - - - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - - 
# - - - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - - 
shade_confidence_intervals   = False
stratify_spectral_components = True  # Change this to False if you don't want to stratify and shade the spectral components
ci = True

# Define the trajbounds for each column of the subplot grid
column_trajbounds = [0, 1]

# Update the components for each row of the subplot grid
row_components = [["U1", "U2", "U3"], ["V1", "V2", "V3"], 
                  ["U1a", "U2a", "U3a"], ["V1a", "V2a", "V3a"],
                  ["Cavgtheta", "S1theta", "S2theta", "wpli_avgtheta"],
                  ["Cavgdelta", "S1delta", "S2delta", "wpli_avgdelta"],
                  ["Cavgripple", "S1ripple", "S2ripple", "wpli_avgripple"]]

                             
component_colors = {
    "U1a": "darkred", "U2a": "red", "U3a": "lightcoral",
    "V1a": "darkblue", "V2a": "blue", "V3a": "lightblue",
    "U1": "darkred", "U2": "red", "U3": "lightcoral",
    "V1": "darkblue", "V2": "blue", "V3": "lightblue",
    "Cavgtheta": "black", "Cavgdelta": "black", "Cavgripple": "black",
    "S1theta": "darkred", "S1delta": "darkred", "S1ripple": "darkred",
    "S2theta": "darkblue", "S2delta": "darkblue", "S2ripple": "darkblue",
    "wpli_avgtheta": "black", "wpli_avgdelta": "black", "wpli_avgripple": "black"
}

component_fill_bases = {
    "U1": 0, "U2": 1, "U3": 2,
    "U1a": 0, "U2a": 1, "U3a": 2,
    "V1": 0, "V2": 1, "V3": 2,
    "V1a": 0, "V2a": 1, "V3a": 2,
    "Cavgtheta": 0, "Cavgdelta": 0, "Cavgripple": 0,
    "S1theta": 1, "S1delta": 1, "S1ripple": 1,
    "S2theta": 1, "S2delta": 1, "S2ripple": 1,
    "wpli_avgtheta": 0, "wpli_avgdelta": 0, "wpli_avgripple": 0
}
# Scale distances
if ci:
    between_scale = 2
else:
    between_scale = 5/8
component_fill_bases = {k: v * between_scale for k, v in component_fill_bases.items()}

scale = 2

line_styles = {
    "U1a": "-", "U2a": "-", "U3a": "-",
    "V1a": "-", "V2a": "-", "V3a": "-",
    "U1": "-", "U2": "-", "U3": "-",
    "V1": "-", "V2": "-", "V3": "-",
    "Cavgtheta": "-", "Cavgdelta": "-", "Cavgripple": "-",
    "S1theta": "-", "S1delta": "-", "S1ripple": "-",
    "S2theta": "-", "S2delta": "-", "S2ripple": "-",
    "wpli_avgtheta": "--", "wpli_avgdelta": "--", "wpli_avgripple": "--"
}

# - - - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - - 
# - - - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - - 
# - - - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - - 
print("Loading data...")
intermediate = "midpattern=true"
folder       = f'/Volumes/MATLAB-Drive/Shared/figures/{intermediate}/tables/'
figfolder    = f'/Volumes/MATLAB-Drive/Shared/figures/{intermediate}/lindist_bootstrap_epoch/zscore'
if not os.path.exists(figfolder):
    os.makedirs(figfolder)
# name   = 'ZT2powerccatime'
name   = 'powerccatime'
append = '_50bin_onlyrcols'
datafile = os.path.join(folder, f'{name}.*')
datafile = glob.glob(datafile)[0]
if os.path.splitext(datafile)[1] == '.csv':
    df = pd.read_csv(datafile)
elif os.path.splitext(datafile)[1] == '.parquet':
    df = pd.read_parquet(datafile)
else:
    raise ValueError(f'Unknown file extension: {os.path.splitext(datafile)[1]}')
print("Done.")
def label_epochs(df, time_col='time', time_threshold=15*60):
    # Calculate time difference
    df.loc[:,'time_diff'] = df[time_col].diff()
    # Initialize the epoch column with 0
    df.loc[:,'epoch'] = 0
    # If time difference is more than threshold, increment epoch
    df.loc[df['time_diff'] > time_threshold, 'epoch'] = 1
    # Calculate the cumulative sum of the epoch column
    df.loc[:, 'epoch'] = df['epoch'].cumsum()
    # Drop the time_diff column as we don't need it anymore
    df.drop(columns=['time_diff'], inplace=True)
    return df
# Test the function on your dataframe
df = df.sort_values(['time','animal']).groupby('animal').apply(label_epochs).reset_index()
if not name.endswith("_epoch"):
    name = name + "_epoch"
df.groupby(['animal','epoch']).time.mean().unstack()

# ----------------------------------------------------
# MAKE SOME VARIABLES
# ----------------------------------------------------
df["Cavgtheta_bin"] = pd.qcut(df["Cavgtheta"], 40, labels=False)
for uv in ["U1", "U2", "U3", "V1", "V2", "V3"]:
    df[f"{uv}a"] = np.abs(df[f"{uv}"])

for i,(u,v) in enumerate(zip(("U1", "U2", "U3"), ("V1", "V2", "V3"))):
    # Measure concordance divided by best possible concordance
    u = df[u]
    v = df[v]
    df[f"R{i+1}"] = (u+v)/np.sqrt(2) / np.sqrt(u **2 + v **2)
    df[f"R{i+1}c"] = df[f"R{i+1}"] / (2*np.mean((u.abs(), v.abs()), axis=0))**2
    df[f"R{i+1}off"] = (u-v)/np.sqrt(2) / np.sqrt(u **2 + v **2)
    df[f"R{i+1}mag"] = np.sqrt(u **2 + v **2)
    df[f"R{i+1}_abs"] = df[f"R{i+1}"].abs()
    df[f"R{i+1}_on_vs_off"] = np.log10(df[f"R{i+1}"] / df[f"R{i+1}off"])
    df[f"R{i+1}_on_vs_mag"] = np.log10(df[f"R{i+1}"] / df[f"R{i+1}mag"])

df["Roverall_on_vs_off"] = np.log10(np.sqrt(np.sum(df[["R1", "R2", "R3"]]**2, axis=1)) / np.sqrt(np.sum(df[["R1off", "R2off", "R3off"]]**2, axis=1)))
df["R_overall_abs"] = np.sqrt(np.sum(df[["R1_abs", "R2_abs", "R3_abs"]]**2, axis=1))
df["R_overall"]     = np.sqrt(np.sum(df[["R1", "R2", "R3"]]**2, axis=1))

df = downcast_dataframe(df)
gc.collect()

# ----------------------------------------------------

# Number of bins and bootstrap samples
n_bins = 45
n_bootstrap_samples = 500
# Update the columns to bootstrap
columns_to_bootstrap = ["U1", "U2", "U3", "V1", "V2", "V3",
                        "U1a", "U2a", "U3a", "V1a", "V2a", "V3a",
                        "Cavgtheta", 
                        "Cavgdelta", "Cavgripple", "S1theta", "S1delta", 
                        "S2theta", "wpli_avgtheta", "S1ripple", "S2ripple",
                        "S2delta", "wpli_avgripple", "wpli_avgdelta", "mua",
                        "vel", "accel",
                        "R1", "R2", "R3", 
                        "R1c", "R2c", "R3c", 
                        "R1off", "R2off", "R3off",
                        "R1mag", "R2mag", "R3mag",
                        "R1_abs", "R2_abs", "R3_abs",
                        "R1_on_vs_off", "R2_on_vs_off", "R3_on_vs_off",
                        "R1_on_vs_mag", "R2_on_vs_mag", "R3_on_vs_mag",
                        "Roverall_on_vs_off", "R_overall_abs", "R_overall"]

# Reinitialize the DataFrame to hold the results
# Bin the lindist data
df["lindist_bin"] = pd.cut(df["lindist"], bins=n_bins)
check_existing = False
if check_existing and os.path.exists(os.path.join(folder, f'{name}_bootstrap{append}.parquet')):
    bootstrap_means_combined = pd.read_parquet(os.path.join(folder, f'{name}_bootstrap{append}.parquet'))
    print("memory usage: ", bootstrap_means_combined.memory_usage().sum() / 1e6, "MB")
else:
    bootstrap_means_combined = []


if append.endswith("onlyrcols"):
    columns_to_bootstrap = [c for c in columns_to_bootstrap if c.startswith("R")]
# bootstrap_means_combined.query('column not in @only_r_columns', inplace=True)

# Rerun the bootstrap for each bin and each column, for each trajectory bound
print("Running bootstrap with {} bins and {} bootstrap samples...".format(n_bins, n_bootstrap_samples))
for trajbound in tqdm([0, 1], desc="trajbound", total=2):
    # Filter the data based on trajbound
    data_trajbound = df.query('trajbound == @trajbound')
    # Loop over the unique combinations of column, trajbound, and lindist_bin
    iters = list(itertools.product(columns_to_bootstrap, 
                                   data_trajbound["lindist_bin"].unique().categories,
                                   data_trajbound["epoch"].unique()))
    bootstrap_means_combined_trajbound = []
    gc.collect()
    for (column, bin_label, epoch) in tqdm(iters, desc="column, lindist_bin", total=len(iters)):
        if check_existing:
            print("Checking for existing entries...")
            existing_entries = bootstrap_means_combined.query('column == @column & lindist_bin == @bin_label & epoch == @epoch').shape[0]
            if existing_entries > 0:
                print("Skipping trajbound {}, column {}, lindist_bin {}, epoch {}".format(trajbound, column, bin_label, epoch))
                continue
        # Get the data for this column, trajbound, and bin_label
        data = data_trajbound.loc[(data_trajbound["lindist_bin"] == bin_label) &
                                  (data_trajbound["epoch"] == epoch),
                                  ["animal",column]]
        # data = data_trajbound.query('lindist_bin == @bin_label & epoch == @epoch')[["animal", column]]
        
        # Find the minimum number of points available for each animal
        min_points_per_animal = data.groupby("animal").size().min()
        
        # Generate bootstrap samples
        bms = []
        for iboot in range(n_bootstrap_samples):
            # Initialize an empty list to hold the data for this bootstrap sample
            bootstrap_sample = []
            
            # Sample the minimum number of points from each animal's data
            for animal, animal_data in data.groupby("animal"):
                sample = (animal_data.sample(min_points_per_animal, replace=True))
                
                # Compute the mean of the bootstrap sample
                bootstrap_mean = sample.mean(numeric_only=True).astype(float)
                bootstrap_var  = sample.var(numeric_only=True).astype(float)
                
                # Add the result to the DataFrame
                bms.append({
                    "iboot": iboot,
                    "lindist_bin": bin_label,
                    "column": column,
                    "bootstrap_mean": bootstrap_mean,
                    "bootstrap_var": bootstrap_var,
                    "trajbound": trajbound,
                    "animal": animal,
                    "epoch": epoch
                })
        bms = pd.DataFrame(bms)
        gc.collect()
        bootstrap_means_combined_trajbound.append(bms)
    bootstrap_means_combined.append(downcast_dataframe(pd.concat(bootstrap_means_combined_trajbound, axis=0)))

# Convert the list of dictionaries to a DataFrame
print("Converting to DataFrame...")
bootstrap_means_combined = pd.concat(bootstrap_means_combined, axis=0)
# Create a new lindist_bin_ind column
print("Creating lindist_bin_ind column...")
bootstrap_means_combined["lindist_bin_ind"] = bootstrap_means_combined["lindist_bin"].apply(lambda x: x.right)
bootstrap_means_combined.loc[:,'bootstrap_mean'] = bootstrap_means_combined.bootstrap_mean.astype(float)
print("Saving to parquet...")
bootstrap_means_combined.to_parquet(os.path.join(folder, f'{name}_bootstrap{append}.parquet'), index=False)
bootstrap_means_combined = pd.read_parquet(os.path.join(folder, f'{name}_bootstrap{append}.parquet'))

# Smooth
print("sorting...")
bootstrap_means_combined.sort_values(by=["animal", "epoch", "column", "trajbound", "iboot", "lindist_bin_ind"], inplace=True)
print("rolling mean...")
bootstrap_means_combined["bootstrap_mean_smooth"] = bootstrap_means_combined.groupby(["animal","epoch","column", "trajbound", "iboot"])["bootstrap_mean"].transform(lambda x: x.rolling(7, 1).mean())
print("rolling interp...")
bootstrap_means_combined["bootstrap_mean_smooth"] = bootstrap_means_combined.groupby(["animal","epoch","column", "trajbound", "iboot"])["bootstrap_mean_smooth"].transform(lambda x: x.interpolate())
# Compute the rolling variance
bootstrap_means_combined["bootstrap_meanvar_smooth"] = bootstrap_means_combined.groupby(
    ["animal","epoch","column", "trajbound", "iboot"]
)["bootstrap_mean"].transform(lambda x: x.rolling(7, 1).var())

# Interpolate the rolling variance
bootstrap_means_combined["bootstrap_meanvar_smooth"] = bootstrap_means_combined.groupby(
    ["animal","epoch","column", "trajbound", "iboot"]
)["bootstrap_meanvar_smooth"].transform(lambda x: x.interpolate())

      
# ----------------------------------------------------
# Normalize the bootstrap_mean values
# ----------------------------------------------------
print("Normalizing bootstrap_mean values...")
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
# scaler = MinMaxScaler()
iters = itertools.product(columns_to_bootstrap, bootstrap_means_combined.animal.unique())
for column, animal in tqdm(iters, desc="feature engineering", total=len(columns_to_bootstrap) * bootstrap_means_combined.animal.nunique()):
    # Get the data for this column and animal
    data = \
            bootstrap_means_combined[(bootstrap_means_combined["column"] == column) & \
            (bootstrap_means_combined["animal"] == animal)]["bootstrap_mean"].values.reshape(-1, 1)
    data_smooth  = \
            bootstrap_means_combined[(bootstrap_means_combined["column"] == column) & \
            (bootstrap_means_combined["animal"] == animal)]["bootstrap_mean_smooth"].values.reshape(-1, 1)
    # Scale the data
    scaled_data        = scaler.fit_transform(data)
    scaled_data_smooth = scaler.fit_transform(data_smooth)
    # Update the DataFrame
    bootstrap_means_combined.loc[(bootstrap_means_combined["column"] == column) & \
            (bootstrap_means_combined["animal"] == animal), "bootstrap_mean"] = scaled_data
    bootstrap_means_combined.loc[(bootstrap_means_combined["column"] == column) & \
            (bootstrap_means_combined["animal"] == animal), "bootstrap_mean_smooth"] = scaled_data_smooth
# bootstrap_means_combined.head()
bootstrap_means_combined.to_parquet(
        os.path.join(folder, f'{name}_bootstrap_normalized{append}_{str(scaler).replace("()","").lower()}.parquet'), index=False)
bootstrap_means_combined = pd.read_parquet(
        os.path.join(folder, f'{name}_bootstrap_normalized{append}_{str(scaler).replace("()","").lower()}.parquet'))

# Add the lindist_bin_mid column back to the DataFrame
print("Adding lindist_bin_mid column...")
bootstrap_means_combined["lindist_bin_mid"] = \
    bootstrap_means_combined["lindist_bin"].apply(lambda x: x.mid)
# FIXME: PROBLEM WITH EPOCH CUTTING, comment out when fix
bootstrap_means_combined.query('animal !="ER1"', inplace=True)
bootstrap_means_combined.to_parquet(
        os.path.join(folder, f'{name}_bootstrap_normalized{append}_{str(scaler).replace("()","").lower()}.parquet'), index=False)


print("Balancing animals...")
animal_bootstrap_means_combined = bootstrap_means_combined.groupby(
        ["iboot", "epoch", "column", "trajbound", "lindist_bin_mid"]).mean().reset_index()
print("Saving to parquet...")
animal_bootstrap_means_combined.to_parquet(
        os.path.join(folder, f'{name}_bootstrap_normalized{append}_balancedanim.parquet'), index=False)

# ----------------------------------------------------
# Read Parquet
# ----------------------------------------------------
bootstrap_means_combined        = pd.read_parquet(os.path.join(folder, f'{name}_bootstrap_normalized{append}.parquet'))
animal_bootstrap_means_combined = pd.read_parquet(os.path.join(folder, f'{name}_bootstrap_normalized{append}_balancedanim.parquet'))

# -------------------------Checking for duplicates------------------- ------ ##

print("Plot the bootstrap means overall")

component_fill_bases = {
    "U1": 0, "U2": 1, "U3": 2,
    "U1a": 0, "U2a": 1, "U3a": 2,
    "V1": 0, "V2": 1, "V3": 2,
    "V1a": 0, "V2a": 1, "V3a": 2,
    "Cavgtheta": 0, "Cavgdelta": 0, "Cavgripple": 0,
    "S1theta": 1, "S1delta": 1, "S1ripple": 1,
    "S2theta": 1, "S2delta": 1, "S2ripple": 1,
    "wpli_avgtheta": 0, "wpli_avgdelta": 0, "wpli_avgripple": 0
}
# Scale distances
scale = 1
skipripple = True
if ci:
    between_scale = 5/8 
else:
    between_scale = 5/8
component_fill_bases = {k: v * between_scale for k, v in component_fill_bases.items()}

rows = 4 if skipripple else 5

for epoch in tqdm(animal_bootstrap_means_combined["epoch"].unique(), desc="epoch"):

    D = animal_bootstrap_means_combined[animal_bootstrap_means_combined["epoch"] == epoch]

    fig, axes = plt.subplots(nrows=rows, ncols=2, figsize=(15, 25))
    smooth = True
    field = "bootstrap_mean_smooth" if smooth else "bootstrap_mean"
    # Set the overall title
    fig.suptitle('Epoch=' + str(epoch), fontsize=20)
    for i, components in enumerate(row_components[:rows]):
        for j, trajbound in enumerate(column_trajbounds):
            # Get the data for this subplot
            data = D[D["column"].isin(components) & (D["trajbound"] == trajbound)]
            print("select components", components, "and unique components", data["column"].unique())
            print("select trajbound", trajbound, "and unique trajbounds", data["trajbound"].unique())
            # Create the subplot
            for component in components:
                if skipripple and "ripple" in component:
                    continue
                component_data = data[data["column"] == component].sort_values(by="lindist_bin_mid")
                # Add the base value to the bootstrap_mean if stratification is active
                # if stratify_spectral_components and component.startswith(("S1", "S2","U", "V")):
                component_data[field] += component_fill_bases[component]
                component_data[field] *= scale
                # Plot the curve
                sns.lineplot(x="lindist_bin_mid", y=field, data=component_data, color=component_colors[component], 
                             errorbar="se", 
                             ax=axes[i, j], linestyle=line_styles[component])
                axes[i,j].axhline(y=scale*(component_fill_bases[component]+0.5), color=component_colors[component], linestyle="dotted")
                if ci:
                    ci_99 = component_data.groupby("lindist_bin_mid").quantile(0.99)[field]
                    ci_1 = component_data.groupby("lindist_bin_mid").quantile(0.01)[field]
                    axes[i,j].plot(ci_99.index, ci_99, color=component_colors[component], linestyle="dashed")
                    axes[i,j].plot(ci_1.index, ci_1, color=component_colors[component], linestyle="dashed")
                if i == 2 or i == 3:
                    axes[i,j].set_ylim([0.3, 1.4])
                # Fill under the curve based on the component type and the flag
                # if component in ["U1", "U2", "U3", "V1", "V2", "V3"] or (stratify_spectral_components and component.startswith(("S1", "S2"))):
                #     axes[i, j].fill_between(component_data["lindist_bin_mid"], 
                #                             component_fill_bases[component], 
                #                             component_data[field], 
                #                             color=component_colors[component], 
                #                             alpha=0.3)
                # Shade the confidence interval around the curve if the flag is set
                # if shade_confidence_intervals:
                #     component_data_bootstrap_samples = component_data[field].tolist()
                #     confidence_interval = np.percentile(component_data_bootstrap_samples, [2.5, 97.5])
                #     axes[i, j].fill_between(component_data["lindist_bin_mid"], confidence_interval[0], confidence_interval[1], color=component_colors[component], alpha=0.1)
            # Set the title and labels
            axes[i, j].set_title(f'Trajbound = {trajbound}')
            axes[i, j].set_xlabel("Lindist Bin Midpoint")
            axes[i, j].set_ylabel("Normalized Bootstrap Mean")
            # Rotate the x-axis labels for readability
            axes[i, j].tick_params(axis='x', rotation=45)
            # Add a legend
            axes[i, j].legend()
    if not os.path.exists(figfolder):
        os.makedirs(figfolder)
    fig.savefig(figfolder + f'lindist_bootstrap{append}_{field}_balancedanim_collapseanimals_{epoch}.png', dpi=300)
    fig.savefig(figfolder + f'lindist_bootstrap{append}_{field}_balancedanim_collapseanimals_{epoch}.svg', dpi=300)
    fig.savefig(figfolder + f'lindist_bootstrap{append}_{field}_balancedanim_collapseanimals_{epoch}.pdf', dpi=300)

# Improve the layout
plt.tight_layout()
# Add space for the overall title
# fig.subplots_adjust(top=0.92)
plt.show()


# ------------------------------------------------------

def plot_with_bootstrap_ci(data, x, y, ax, color, label):
    # Extract unique values of x and sort them
    x_values = sorted(data[x].unique())
    # For each unique x value, compute the mean and confidence interval of y
    means = []
    lower_bound = []
    upper_bound = []
    for x_val in x_values:
        y_samples = data[data[x] == x_val][y]
        means.append(np.mean(y_samples))
        lower_bound.append(np.percentile(y_samples, 2.5))
        upper_bound.append(np.percentile(y_samples, 97.5))
    # Plot the mean
    ax.plot(x_values, means, color=color, label=label)
    # Shade the confidence interval
    ax.fill_between(x_values, lower_bound, upper_bound, color=color, alpha=0.3)

# Example usage
# fig, ax = plt.subplots()
# plot_with_bootstrap_ci(group.query('trajbound == 0'), "lindist_bin_mid", "your_y_field_here", ax, "blue", "label_here")
# plt.show()


def plot_by_epoch(df, column, epochs, field="bootstrap_mean_smooth"):
    # Get the unique epochs
    unique_epochs = sorted(df[epochs].unique())
    
    # Create colormap based on the number of epochs
    colors = cm.coolwarm(np.linspace(0, 1, len(unique_epochs)))
    
    # Create a dictionary to map epoch to color
    epoch_color_map = dict(zip(unique_epochs, colors))
    
    # Create a subplot
    fig, ax = plt.subplots(2,1,figsize=(10, 7))
    
    # Group data by epoch and plot each group
    for epoch, group in df.groupby(epochs):
        # Subset to the column of interest
        group = group[group["column"] == column]
        # Sort values for consistent plotting
        group = group.sort_values(by="lindist_bin_mid")
        # Plot the line for this epoch
        # sns.lineplot(x="lindist_bin_mid", y=field, 
        #              data=group.query('trajbound == 0'), 
        #              ax=ax[0],
        #              color=epoch_color_map[epoch], label=f'Epoch {epoch}')
        plot_with_bootstrap_ci(group.query('trajbound == 0'),
                               "lindist_bin_mid", field, ax[0],
                               epoch_color_map[epoch], f'E{epoch}')
        # sns.lineplot(x="lindist_bin_mid", y=field, 
        #              data=group.query('trajbound == 1'), 
        #              ax=ax[1],
        #              color=epoch_color_map[epoch], label=f'Epoch {epoch}')
        plot_with_bootstrap_ci(group.query('trajbound == 1'),
                               "lindist_bin_mid",field, ax[1],
                               epoch_color_map[epoch], f'')
    
    # Set title and labels
    fig.suptitle(f'{column} over epochs')
    fig.legend()
    ax[1].set_xlabel('lindist_bin_mid')
    ax[1].set_ylabel(field)
    ax[0].set_ylabel(field)
    
    # Show the plot
    plt.show()

print("Plotting by epoch")
# Call the function
for column in tqdm(component_fill_bases.keys(),total=len(component_fill_bases.keys())):
    plot_by_epoch(animal_bootstrap_means_combined, column, "epoch")
    plt.savefig(os.path.join(figfolder, f"{column}_bootstrap_mean_smooth_by_epoch.png"), dpi=300)
    plt.savefig(os.path.join(figfolder, f"{column}_bootstrap_mean_smooth_by_epoch.pdf"), dpi=300)
    # plt.close('all')

# ------------------------------------------------------
# By animal
# ------------------------------------------------------
def plot_by_animal(df, column, field="bootstrap_mean_smooth"):
    # Get the unique animals, epochs, and trajbounds
    unique_animals = sorted(df["animal"].unique())
    unique_epochs = sorted(df["epoch"].unique())
    unique_trajbounds = sorted(df["trajbound"].unique())
    
    # Create colormap based on the number of epochs
    colors = cm.coolwarm(np.linspace(0, 1, len(unique_epochs)))

    # Create a dictionary to map epoch to color
    epoch_color_map = dict(zip(unique_epochs, colors))
    
    # Create a subplot
    fig, axes = plt.subplots(len(unique_animals), len(unique_trajbounds),
                             figsize=(10, 7 * len(unique_animals)))

    # Group data by animal, epoch, and trajbound and plot each group
    for animal in tqdm(unique_animals, desc="Animals", total=len(unique_animals)):
        for trajbound in unique_trajbounds:
            for epoch in unique_epochs:
                # Subset to the data of interest
                data = df[(df["animal"] == animal) & (df["trajbound"] == trajbound) & (df["epoch"] == epoch) & (df["column"] == column)]
                # Sort values for consistent plotting
                data = data.sort_values(by="lindist_bin_mid")
                # Plot the line for this animal, epoch, and trajbound
                if not data.empty:
                    # sns.lineplot(
                    #     x="lindist_bin_mid", y=field, data=data,
                    #     ax=axes[unique_animals.index(animal), unique_trajbounds.index(trajbound)],
                    #     color=epoch_color_map[epoch], label=f'Epoch {epoch}'
                    # )
                    plot_with_bootstrap_ci(data, "lindist_bin_mid", field,
                                           ax=axes[unique_animals.index(animal),
                                                   unique_trajbounds.index(trajbound)],
                                           color=epoch_color_map[epoch],
                                           label=f'E{epoch}' if trajbound == 0 else '')
                    # Set title and labels
                    axes[unique_animals.index(animal), unique_trajbounds.index(trajbound)].set_title(f'Animal {animal}, Trajbound {trajbound}')
                    axes[unique_animals.index(animal), unique_trajbounds.index(trajbound)].set_xlabel('lindist_bin_mid')
                    axes[unique_animals.index(animal), unique_trajbounds.index(trajbound)].set_ylabel(field)

    # Show the plot
    plt.tight_layout()
    plt.show()

print("Plotting by animal")
# The test function call has been commented out because we don't currently have access to the bootstrap_means_combined DataFrame
for column in tqdm(component_fill_bases.keys(),total=len(component_fill_bases.keys())):
    plot_by_animal(bootstrap_means_combined, column)
    plt.savefig(os.path.join(figfolder, f"{column}_bootstrap_mean_smooth_by_animal.png"), dpi=300)
    plt.savefig(os.path.join(figfolder, f"{column}_bootstrap_mean_smooth_by_animal.pdf"), dpi=300)
    # plt.close('all')


# ------------------------------------------------------
# RAW
# ------------------------------------------------------

def plot_property_vs_lindist(df, property_name, N=100_000):
    """
    Plot the given property against lindist for both trajbound values.
    
    Parameters:
    - df: DataFrame containing the data
    - property_name: Name of the property/column to plot
    - N: Number of points to subsample
    """
    
    # Sort the data by lindist
    
    # Set up the figure and axis
    fig, ax = plt.subplots(2,1,figsize=(12, 6))
    ax = ax.flatten()
    
    # Plot for each trajbound
    for i, trajbound in enumerate([0, 1]):
        traj_data = df[df["trajbound"] == trajbound]
        
        # Subsample
        sampled_data = traj_data.sample(n=min(N, len(traj_data)), random_state=42)
        sampled_data = sampled_data.sort_values(by=["trajbound","lindist"])
        
        # Plot
        ax[i].axhline(y=0, color="red", linestyle="dashed", alpha=0.2)
        ax[i].axvline(x=0.4, color="black", linestyle="dashed", label="Turn 1")
        ax[i].axvline(x=0.6, color="black", linestyle="dashed", label="Turn 2")
        ax[i].plot(sampled_data["lindist"], sampled_data[property_name], 
                 linewidth=0.2, alpha=0.5)
        ax[i].set_title(f"Trajbound {trajbound}")
    
        # Finalize the plot
        plt.xlabel("lindist")
        plt.ylabel(property_name)
        plt.legend()
        plt.grid(True)
    fig.suptitle(f"{property_name} vs lindist")
    fig.savefig(os.path.join(figfolder, f"{property_name}_vs_lindist.png"), dpi=300)
    plt.show()

plot_property_vs_lindist(df, "U1")
plot_property_vs_lindist(df, "U2")
plot_property_vs_lindist(df, "U3")

plot_property_vs_lindist(df, "V1")
plot_property_vs_lindist(df, "V2")
plot_property_vs_lindist(df, "V3")

plot_property_vs_lindist(df, "R1")
plot_property_vs_lindist(df, "R2")
plot_property_vs_lindist(df, "R3")

plot_property_vs_lindist(df, "Cavgtheta")

import os
import matplotlib.pyplot as plt

def plot_property_vs_lindist_by_animal(df, property_name, N=100_000, winsize=None, method="mean"):
    """
    Plot the given property against lindist for both trajbound values and for each animal.
    
    Parameters:
    - df: DataFrame containing the data
    - property_name: Name of the property/column to plot
    - N: Number of points to subsample
    - winsize: Window size for rolling mean. If None, no rolling mean is applied.
    - figfolder: Folder to save the figures
    """
    
    # Set of unique animals
    animals = df['animal'].unique()
    
    # Set up the figure and axes
    fig, axs = plt.subplots(len(animals), 2, figsize=(12, 6 * len(animals)))
    
    # Plot for each trajbound and each animal
    for i, animal in enumerate(animals):
        animal_data = df[df["animal"] == animal]
        
        for j, trajbound in enumerate([0, 1]):
            traj_data = animal_data[animal_data["trajbound"] == trajbound]
            
            # Subsample
            sampled_data = traj_data.sample(n=min(N, len(traj_data)), random_state=42)
            sampled_data = sampled_data.sort_values(by=["trajbound", "lindist"])
            
            if winsize:
                if method == "mean":
                    sampled_data[property_name] = sampled_data[property_name].rolling(winsize, min_periods=1).mean().interpolate()
                elif method == "var":
                    print("variance")
                    sampled_data[property_name] = sampled_data[property_name].rolling(winsize, min_periods=1).var().rolling(7).mean().interpolate()
                elif method == "std":
                    print("std")
                    sampled_data[property_name] = sampled_data[property_name].rolling(winsize, min_periods=1).std().interpolate()
                elif method == "median":
                    print("median")
                    sampled_data[property_name] = sampled_data[property_name].rolling(winsize, min_periods=1).median().interpolate()
        
            # Plot
            axs[i][j].axhline(y=0, color="red", linestyle="dashed", alpha=0.2)
            axs[i][j].axvline(x=0.4, color="black", linestyle="dashed", label="Turn 1")
            axs[i][j].axvline(x=0.6, color="black", linestyle="dashed", label="Turn 2")
            axs[i][j].plot(sampled_data["lindist"], sampled_data[property_name], linewidth=0.2 if winsize is None else 2, alpha=0.5 if winsize is None else None)
            axs[i][j].set_title(f"Animal {animal} | Trajbound {trajbound}")
            
            # Labels and legend
            axs[i][j].set_xlabel("lindist")
            axs[i][j].set_ylabel(property_name)
            axs[i][j].legend()
            axs[i][j].grid(True)
    
    # Overall title
    if not method.startswith("\n"):
        method = "\n" + method
    fig.suptitle(f"{property_name} vs lindist by Animal{method}", y=1.02)
    fig.subplots_adjust(top=0.95)
    
    # Save the figure
    figname = f"{property_name}_vs_lindist_by_animal"
    if winsize:
        figname += f"_win={winsize}"
    if "mean" in method:
        method = ""
    else:
        method = "_" + method.replace("\n", "")
    print("Saving " + figname + method)
    fig.savefig(os.path.join(figfolder, f"{figname}{method}.png"), dpi=300)
    fig.savefig(os.path.join(figfolder, f"{figname}{method}.pdf"), dpi=300)
    plt.tight_layout()
    fig.subplots_adjust(top=0.95)
    plt.show()


plot_property_vs_lindist_by_animal(df, "U1")
plot_property_vs_lindist_by_animal(df, "U2")
plot_property_vs_lindist_by_animal(df, "U3")

plot_property_vs_lindist_by_animal(df, "V1")
plot_property_vs_lindist_by_animal(df, "V2")
plot_property_vs_lindist_by_animal(df, "V3")

plot_property_vs_lindist_by_animal(df, "R1")
plot_property_vs_lindist_by_animal(df, "R2")
plot_property_vs_lindist_by_animal(df, "R3")

plot_property_vs_lindist_by_animal(df, "Cavgtheta")
plot_property_vs_lindist_by_animal(df, "Cavgdelta")
plot_property_vs_lindist_by_animal(df, "vel")
plot_property_vs_lindist_by_animal(df, "accel")

plot_property_vs_lindist_by_animal(df, "R1", winsize=100)
plot_property_vs_lindist_by_animal(df, "R1", winsize=100, method="std")

plot_property_vs_lindist_by_animal(df, "R1d")
plot_property_vs_lindist_by_animal(df, "R1d", winsize=100)
plot_property_vs_lindist_by_animal(df, "R1d", winsize=100, method="std")
plot_property_vs_lindist_by_animal(df, "R2d")
plot_property_vs_lindist_by_animal(df, "R2d", winsize=100)
plot_property_vs_lindist_by_animal(df, "R2d", winsize=100, method="std")
plot_property_vs_lindist_by_animal(df, "R3d")
plot_property_vs_lindist_by_animal(df, "R3d", winsize=100)


pg = sns.PairGrid(df.sample(10_000), vars=["R1c", "R2c", "R3c"], hue="trajbound", palette="Set2")
pg.map_diag(sns.histplot)
pg.map_upper(sns.scatterplot)
# pg.map_lower(sns.hexbinplot)
pg.add_legend()

pg = sns.PairGrid(df.sample(10_000), vars=["U1", "U2", "U3", "V1", "V2", "V3"], hue="trajbound", palette="Set2")
pg.map_diag(sns.histplot)
pg.map_upper(sns.scatterplot)

# WARNING: TODO we should qqplot the cavgtheta and r1/2/3
# and do so by epoch as well

def assign_lindistzone(ld):
    if ld < 0.4:
        return "center"
    elif ld < 0.6:
        return "turn"
    else:
        return "side"
df['lindistzone'] = df['lindist'].apply(assign_lindistzone)


g=sns.relplot(data=df.sample(10_000), x="epoch", y="R1d", row="trajbound",
           palette="coolwarm", kind="line", err_style="bars", ci=95)
g.savefig(os.path.join(figfolder, f"R1d_by_epoch.png"), dpi=300)
g=sns.relplot(data=df.sample(10_000), x="epoch", y="R2d", row="trajbound",
            palette="coolwarm", kind="line", err_style="bars", ci=95)
g.savefig(os.path.join(figfolder, f"R2d_by_epoch.png"), dpi=300)
g=sns.relplot(data=df.sample(10_000), x="epoch", y="R3d", row="trajbound",
            palette="coolwarm", kind="line", err_style="bars", ci=95)
g.savefig(os.path.join(figfolder, f"R3d_by_epoch.png"), dpi=300)


g=sns.relplot(data=df.sample(10_000), x="epoch", y="R1d", row="trajbound", hue="lindistzone", hue_order=["center","turn","side"],
           palette="coolwarm", kind="line", err_style="bars", ci=95)
g.map(lambda **kwargs: plt.axhline(0, color="black", linestyle="dashed"))
g.savefig(os.path.join(figfolder, f"R1d_by_epoch_by_zone.png"), dpi=300)
g.savefig(os.path.join(figfolder, f"R1d_by_epoch_by_zone.pdf"), dpi=300)
g=sns.relplot(data=df.sample(10_000), x="epoch", y="R2d", row="trajbound", hue="lindistzone", hue_order=["center","turn","side"],
            palette="coolwarm", kind="line", err_style="bars", ci=95)
g.map(lambda **kwargs: plt.axhline(0, color="black", linestyle="dashed"))
g.savefig(os.path.join(figfolder, f"R2d_by_epoch_by_zone.png"), dpi=300)
g.savefig(os.path.join(figfolder, f"R2d_by_epoch_by_zone.pdf"), dpi=300)
g=sns.relplot(data=df.sample(10_000), x="epoch", y="R3d", row="trajbound",hue="lindistzone", hue_order=["center","turn","side"],
            palette="coolwarm", kind="line", err_style="bars", ci=95)
g.map(lambda **kwargs: plt.axhline(0, color="black", linestyle="dashed"))
g.savefig(os.path.join(figfolder, f"R3d_by_epoch_by_zone.png"), dpi=300)
g.savefig(os.path.join(figfolder, f"R3d_by_epoch_by_zone.pdf"), dpi=300)

# ------------------------------------------------------

g=sns.relplot(data=df.sample(80_000), x="epoch", y="R1d", row="trajbound", hue="lindistzone", hue_order=["center","turn","side"],
              col="animal", 
           palette="coolwarm", kind="line", err_style="bars", ci=95)
g.map(lambda **kwargs: plt.axhline(0, color="black", linestyle="dashed"))
g.savefig(os.path.join(figfolder, f"R1d_by_epoch_by_zone_by_animal.png"), dpi=300)
g.savefig(os.path.join(figfolder, f"R1d_by_epoch_by_zone_by_animal.pdf"), dpi=300)
g=sns.relplot(data=df.sample(10_000), x="epoch", y="R2d", row="trajbound", hue="lindistzone", hue_order=["center","turn","side"],
              col="animal",
            palette="coolwarm", kind="line", err_style="bars", ci=95)
g.map(lambda **kwargs: plt.axhline(0, color="black", linestyle="dashed"))
g.savefig(os.path.join(figfolder, f"R2d_by_epoch_by_zone_by_animal.png"), dpi=300)
g.savefig(os.path.join(figfolder, f"R2d_by_epoch_by_zone_by_animal.pdf"), dpi=300)
g=sns.relplot(data=df.sample(10_000), x="epoch", y="R3d", row="trajbound",hue="lindistzone", hue_order=["center","turn","side"],
              col="animal",
            palette="coolwarm", kind="line", err_style="bars", ci=95)
g.map(lambda **kwargs: plt.axhline(0, color="black", linestyle="dashed"))
g.savefig(os.path.join(figfolder, f"R3d_by_epoch_by_zone_by_animal.png"), dpi=300)
g.savefig(os.path.join(figfolder, f"R3d_by_epoch_by_zone_by_animal.pdf"), dpi=300)

# ------------------------------------------------------

df["R1d_abs"] = df["R1d"].abs()
g=sns.relplot(data=df.sample(80_000), 
              x="epoch", y="R1d_abs",
              row="trajbound", hue="lindistzone",
              hue_order=["center","turn","side"], col="animal",
              palette="coolwarm", kind="line", err_style="bars",
              ci=95)
g.map(lambda **kwargs: plt.axhline(0, color="black",
                                   linestyle="dashed"))

# ------------------------------------------------------
# Create 10 quantile bins of Cavgtheta
# ------------------------------------------------------



for col in ["R1d_abs", "R2d_abs", "R3d_abs", "R_overall_abs"]:
    # Relation of Cavgtheta to R1d
    g=sns.relplot(data=df,
                  x="Cavgtheta_bin", y=col,
                  row="trajbound", hue="lindistzone",
                  hue_order=["center","turn","side"], 
                  palette="coolwarm", kind="line", err_style="bars",
                  ci=95)
    g.fig.suptitle(f"{col} vs Cavgtheta_bin mean")
    g.fig.subplots_adjust(top=0.95)
    g.fig.savefig(os.path.join(figfolder, f"{col}_vs_Cavgtheta_bin_mean.png"))
    g.fig.savefig(os.path.join(figfolder, f"{col}_vs_Cavgtheta_bin_mean.pdf"))
    # Relation of Cavgtheta to R1d
    g=sns.relplot(data=df,
                  x="Cavgtheta_bin", y=col,
                  row="trajbound", hue="lindistzone",
                  hue_order=["center","turn","side"], 
                  estimator="std",
                  palette="coolwarm", kind="line", err_style="bars",
                  ci=95)
    g.fig.suptitle(f"{col} vs Cavgtheta_bin std")
    g.fig.subplots_adjust(top=0.95)
    g.fig.savefig(os.path.join(figfolder, f"{col}_vs_Cavgtheta_bin_std.png"))
    g.fig.savefig(os.path.join(figfolder, f"{col}_vs_Cavgtheta_bin_std.pdf"))


g=sns.pairplot(df.sample(5_000), vars=["R1d_abs", "R2d_abs", "R3d_abs", "R_overall_abs", "Cavgtheta", "wpli_avgtheta"], hue="trajbound", palette="Set2")
g.map_diag(sns.histplot)
# sns.map_upper(sns.scatterplot, alpha=0.1)
g.map_lower(sns.regplot, scatter_kws={"alpha":0.1})
plt.savefig(os.path.join(figfolder, f"PairPlot_R1dR2dR3dRoverall_vs_Cavgtheta.png"), dpi=300)
plt.savefig(os.path.join(figfolder, f"PairPlot_R1dR2dR3dRoverall_vs_Cavgtheta.pdf"), dpi=300)

g=sns.pairplot(df.sample(5_000), vars=["R1d_abs", "R2d_abs", "R3d_abs", "R_overall_abs", "Cavgtheta", "wpli_avgtheta"], hue="epoch", palette="coolwarm", kind=None, diag_kind=None)
g.map_diag(sns.histplot)
g.map_upper(sns.kdeplot)
g.map_lower(sns.regplot, scatter_kws={"alpha":0.1})
plt.savefig(os.path.join(figfolder, f"PairPlot_R1dR2dR3dRoverall_vs_Cavgtheta_epoch.png"), dpi=300)
plt.savefig(os.path.join(figfolder, f"PairPlot_R1dR2dR3dRoverall_vs_Cavgtheta_epoch.pdf"), dpi=300)

g=sns.pairplot(df.sample(5_000), vars=["R1d", "R2d", "R3d", "R_overall", "Cavgtheta", "wpli_avgtheta","S1theta","S2theta"], hue="epoch", palette="coolwarm", kind=None, diag_kind=None)
g.map_diag(sns.histplot)
g.map_upper(sns.kdeplot)
g.map_lower(sns.regplot, scatter_kws={"alpha":0.1})
plt.savefig(os.path.join(figfolder, f"PairPlot_nonabs_R1dR2dR3dRoverall_vs_Cavgtheta_epoch.png"), dpi=300)
plt.savefig(os.path.join(figfolder, f"PairPlot_nonabs_R1dR2dR3dRoverall_vs_Cavgtheta_epoch.pdf"), dpi=300)


# TODO: ON VS OFF
tmp=df.sample(10_000).query("Roverall_on_vs_off < 8 & R1_on_vs_off < 8 & R2_on_vs_off < 8 & R3_on_vs_off < 8")
g=sns.pairplot(tmp, vars=["R1_on_vs_off", "R2_on_vs_off", "R3_on_vs_off", "Roverall_on_vs_off", "Cavgtheta", "wpli_avgtheta"], hue="epoch", palette="coolwarm", kind=None, diag_kind=None)
g.map_diag(sns.histplot)
g.map_lower(sns.regplot, scatter_kws={"alpha":0.1})
g.map_upper(sns.kdeplot)


# Let's catplot(kind='bar') the level of on_vs_off for each epoch (x) and trajbound (hue) and animal (row)

tmp = df.sample(100_000).query("Roverall_on_vs_off < 8 & R1_on_vs_off < 8 & R2_on_vs_off < 8 & R3_on_vs_off < 8")
g=sns.catplot(data=tmp, x="epoch", y="Roverall_on_vs_off", hue="trajbound", col="animal", kind="bar", palette="Set2")
g.savefig(os.path.join(figfolder, f"Roverall_on_vs_off_by_epoch_by_animal.png"), dpi=300)
g.savefig(os.path.join(figfolder, f"Roverall_on_vs_off_by_epoch_by_animal.pdf"))
g=sns.catplot(data=tmp, x="epoch", y="R1_on_vs_off", hue="trajbound", col="animal", kind="bar", palette="Set2")
g.savefig(os.path.join(figfolder, f"R1_on_vs_off_by_epoch_by_animal.png"), dpi=300)
g.savefig(os.path.join(figfolder, f"R1_on_vs_off_by_epoch_by_animal.pdf"))
g=sns.catplot(data=tmp, x="epoch", y="R2_on_vs_off", hue="trajbound", col="animal", kind="bar", palette="Set2")
g.savefig(os.path.join(figfolder, f"R2_on_vs_off_by_epoch_by_animal.png"), dpi=300)
g.savefig(os.path.join(figfolder, f"R2_on_vs_off_by_epoch_by_animal.pdf"))
g=sns.catplot(data=tmp, x="epoch", y="R3_on_vs_off", hue="trajbound", col="animal", kind="bar", palette="coolwarm")
g.savefig(os.path.join(figfolder, f"R3_on_vs_off_by_epoch_by_animal.png"), dpi=300)
g.savefig(os.path.join(figfolder, f"R3_on_vs_off_by_epoch_by_animal.pdf"))
