import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import statsmodels.api as sm
import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
from sklearn.utils import resample
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from tqdm import tqdm
import itertools
import matplotlib.cm as cm

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
folder = f'/Volumes/MATLAB-Drive/Shared/figures/{intermediate}/tables/'
figfolder = f'/Volumes/MATLAB-Drive/Shared/figures/{intermediate}/ccatime_regress/'
if not os.path.exists(figfolder):
    os.makedirs(figfolder)
# name   = 'ZT2powerccatime'
name   = 'powerccatime'
append = '_50bin'
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

for uv in ["U1", "U2", "U3", "V1", "V2", "V3"]:
    df[f"{uv}a"] = np.abs(df[f"{uv}"])


# List of variables to inspect
variables_to_inspect = ["U1", "U2", "U3", "V1", "V2", "V3", "U1a", "U2a", "U3a", "V1a", "V2a",
                        "V3a", "Cavgtheta", "Cavgdelta", "Cavgripple", "S1theta", "S1delta", "S1ripple", 
                        "S2theta", "S2delta", "S2ripple", "wpli_avgtheta", "wpli_avgdelta", "wpli_avgripple"]

# Results container
regression_results = {}
regression_models = {}

# Scaling method flag: 'standard', 'minmax', or None
scale_method = 'standard'  # change this value as needed

# Quantile filter: set a value between 0 and 0.5, or set to None for no filtering
quantile_filter = 0.01  # change this value as needed

if scale_method == 'standard':
    scaler = StandardScaler()
elif scale_method == 'minmax':
    scaler = MinMaxScaler()
else:
    scaler = None

# Loop over each variable in the inspection list
tmp = df.dropna()
for variable in variables_to_inspect:

    X = tmp[[variable]].values
    y = tmp['mua'].values
     # If quantile filtering is defined, filter the data based on quantiles
    if quantile_filter:
        X_low, X_high = np.quantile(X, [quantile_filter, 1-quantile_filter])
        y_low, y_high = np.quantile(y, [quantile_filter, 1-quantile_filter])
        mask = (X > X_low).ravel() & (X < X_high).ravel() & (y > y_low).ravel() & (y < y_high).ravel()
        X = X[mask]
        y = y[mask]
    if X.shape[0] == 0:
        continue

    # If a scaler is defined, fit and transform the data
    if scaler:
        X = scaler.fit_transform(X)
        y = scaler.fit_transform(y.reshape(-1, 1)).ravel()

    # Add an intercept to the independent variable
    X = sm.add_constant(X)
    
    # Perform the linear regression
    model = sm.OLS(y, X).fit()
    # Store the model summary in the results container
    regression_results[variable] = model.summary()
    regression_models[variable] = model


# You can then access and print the summary for each variable like this:
for variable, summary in regression_results.items():
    print(f"Regression results for {variable}:\n")
    print(summary)
    print("\n" + "="*80 + "\n")



# Initialize the figure
fig, axes = plt.subplots(nrows=6, ncols=4, figsize=(20, 30))
# Flatten the axes for easier indexing
axes = axes.flatten()
# For each variable in the inspection list, perform regression and plot
for idx, variable in enumerate(variables_to_inspect):
    print("Plotting for", variable)
    # Define the independent variable (with an intercept)
    qhigh = df[variable].quantile(0.99)
    qlow = df[variable].quantile(0.01)
    dfsamp  = df.query(f"{variable} < {qhigh} and {variable} > {qlow}").sample(7_500)
    X = sm.add_constant(dfsamp[variable])
    # Define the dependent variable
    y = dfsamp['mua']
    # Perform the linear regression
    model = regression_models[variable]
    
    # Plot the data and the regression line
    # sns.regplot(x=variable, y='mua', data=df, ax=axes[idx], scatter_kws={'s': 0.5}, line_kws={"color": "red"})
    sns.scatterplot(x=variable, y='mua', data=dfsamp, ax=axes[idx], s=0.5)
    sns.lineplot(x=dfsamp[variable], y=model.predict(X), ax=axes[idx], color="red")
    
    # Set title with R^2 and p-value
    title_text = f"{variable} vs. MUA\n" + r"$R^2 = $" + f"{model.rsquared:.2f}, p = {model.f_pvalue:.2e}"
    axes[idx].set_title(title_text)
# Adjust the layout
plt.show()
plt.savefig(os.path.join(figfolder, f"{name}_regress.png"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(figfolder, f"{name}_regress.pdf"), dpi=300, bbox_inches="tight")

# -----------------------------
# EPOCH WISE TOTALS
# -----------------------------

# List of variables to inspect (as provided earlier)
plt.close('all')
variables_to_inspect = ["U1", "U2", "U3", "V1", "V2", "V3", "U1a", "U2a", "U3a", "V1a", "V2a",
                        "V3a", "Cavgtheta", "Cavgdelta", "Cavgripple", "S1theta", "S1delta", "S1ripple", 
                        "S2theta", "S2delta", "S2ripple", "wpli_avgtheta", "wpli_avgdelta", "wpli_avgripple", "mua", "vel", "accel"]

# Create a copy of the original dataframe
df_filtered = df.copy()
# For each variable, filter based on quantiles and then sample 10,000 points
df_filts = []
stand = StandardScaler()
minmax = MinMaxScaler()
for variable in tqdm(variables_to_inspect, desc="Filtering"):
    low, high = df[variable].quantile([0.01, 0.99])
    df_filtered = df_filtered[(df_filtered[variable] > low) & (df_filtered[variable] < high)]
    df_filtered = df_filtered.sample(min(len(df_filtered), 10000), random_state=42)
    tmp =  df_filtered.loc[:, [variable, 'epoch','animal']] 
    tmp.rename(columns={variable: 'value'}, inplace=True)
    tmp.loc[:, 'variable'] = variable
    tmp.loc[:, 'mmvalue'] = tmp.groupby('animal').value.transform(lambda x: minmax.fit_transform(x.values.reshape(-1,1)).ravel())
    tmp.loc[:, 'stvalue'] = tmp.groupby('animal').value.transform(lambda x: stand.fit_transform(x.values.reshape((-1,1))).ravel())
    tmp.query(f"epoch <= 7", inplace=True)
    df_filts.append(tmp)
df_filtered = pd.concat(df_filts)

# Create the catplot
g = sns.catplot(data=df_filtered, x='epoch', y='value', col='variable',
                kind='point', col_wrap=4, color='black', height=3, aspect=1.5,
                sharey=False, col_order=variables_to_inspect)
# Adjust the y-axis labels to show the variable name
g.set_titles("{col_name}")
g.fig.suptitle("Epoch-wise totals")
plt.show()
plt.savefig(os.path.join(figfolder, f"{name}_epochwise.png"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(figfolder, f"{name}_epochwise.pdf"), dpi=300, bbox_inches="tight")

# Create the catplot
g = sns.catplot(data=df_filtered, x='epoch', y='stvalue', col='variable',
                kind='point', col_wrap=4, color='black', height=3, aspect=1.5,
                sharey=False, col_order=variables_to_inspect)
# Adjust the y-axis labels to show the variable name
g.set_titles("{col_name} (standardized)")
g.fig.suptitle("Epoch-wise totals (standardized)")
g.fig.subplots_adjust(top=0.95)
plt.show()
plt.savefig(os.path.join(figfolder, f"{name}_epochwise_standardized.png"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(figfolder, f"{name}_epochwise_standardized.pdf"), dpi=300, bbox_inches="tight")


# Create the catplot
g = sns.catplot(data=df_filtered, x='epoch', y='mmvalue', col='variable',
                kind='point', col_wrap=4, color='black', height=3, aspect=1.5,
                sharey=False, col_order=variables_to_inspect)
# Adjust the y-axis labels to show the variable name
g.set_titles("{col_name} (standardized)")
g.fig.suptitle("Epoch-wise totals (standardized)")
g.fig.subplots_adjust(top=0.95)
plt.show()
plt.savefig(os.path.join(figfolder, f"{name}_epochwise_standardized.png"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(figfolder, f"{name}_epochwise_standardized.pdf"), dpi=300, bbox_inches="tight")

# -----------------------------
# Animal wise version of last three plots (hue)
# -----------------------------

# Create the catplot
g = sns.catplot(data=df_filtered, x='epoch', y='value', col='variable',
                kind='point', col_wrap=4, hue='animal', height=3, aspect=1.5,
                sharey=False, col_order=variables_to_inspect)
# Adjust the y-axis labels to show the variable name
g.set_titles("{col_name}")
g.fig.suptitle("Epoch-wise split by animal")
plt.show()
plt.savefig(os.path.join(figfolder, f"{name}_epochwise_animal.png"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(figfolder, f"{name}_epochwise_animal.pdf"), dpi=300, bbox_inches="tight")

# Create the catplot
g = sns.catplot(data=df_filtered, x='epoch', y='mmvalue', col='variable',
                kind='point', col_wrap=4, height=3, aspect=1.5,
                sharey=False, col_order=variables_to_inspect, hue='animal')
# Adjust the y-axis labels to show the variable name
g.set_titles("{col_name} (minmax)")
g.fig.suptitle("Epoch-wise totals (minmax)")
g.fig.subplots_adjust(top=0.95)
plt.show()
plt.savefig(os.path.join(figfolder, f"{name}_epochwise_minmax_animal.png"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(figfolder, f"{name}_epochwise_standardized_animal.pdf"), dpi=300, bbox_inches="tight")



# Create the catplot
g = sns.catplot(data=df_filtered, x='epoch', y='stvalue', col='variable',
                kind='point', col_wrap=4, height=3, aspect=1.5,
                sharey=False, col_order=variables_to_inspect, hue='animal')
# Adjust the y-axis labels to show the variable name
g.set_titles("{col_name} (standardized)")
g.fig.suptitle("Epoch-wise totals (standardized)")
g.fig.subplots_adjust(top=0.95)
plt.show()
plt.savefig(os.path.join(figfolder, f"{name}_epochwise_standardized_animal.png"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(figfolder, f"{name}_epochwise_standardized_animal.pdf"), dpi=300, bbox_inches="tight")

# -----------------------------
# MUA Controls
# -----------------------------
for variable in variables_to_inspect:
    df[f"{variable}" + "_divmua"] = df[variable] / df['mua']
    df[f"{variable}" + "_divmua"].replace([np.inf, -np.inf], np.nan, inplace=True)

# Create a copy of the original dataframe
df_filtered = df.copy().dropna()
# For each variable, filter based on quantiles and then sample 10,000 points
df_filts = []
stand = StandardScaler()
minmax = MinMaxScaler()
mua_vars = [f"{variable}" + "_divmua" for variable in variables_to_inspect]
for variable in tqdm(mua_vars, desc="Filtering"):
    tmp =  df_filtered.loc[:, [variable, 'epoch','animal']]
    low, high = tmp[variable].quantile([0.01, 0.99])
    tmp = tmp[(tmp[variable] > low) & (tmp[variable] < high)]
    tmp = tmp.sample(min(len(tmp), 10000), random_state=42)
    if tmp.empty:
        continue
    tmp.rename(columns={variable: 'value'}, inplace=True)
    tmp.loc[:, 'variable'] = variable
    tmp.loc[:, 'mmvalue'] = tmp.groupby('animal').value.transform(lambda x: minmax.fit_transform(x.values.reshape(-1,1)).ravel())
    tmp.loc[:, 'stvalue'] = tmp.groupby('animal').value.transform(lambda x: stand.fit_transform(x.values.reshape((-1,1))).ravel())
    tmp.query(f"epoch <= 7", inplace=True)
    df_filts.append(tmp)
df_filtered = pd.concat(df_filts)

# Create the catplot
g = sns.catplot(data=df_filtered, x='epoch', y='value', col='variable',
                kind='point', col_wrap=4, color='black', height=3, aspect=1.5,
                sharey=False, col_order=variables_to_inspect)
# Adjust the y-axis labels to show the variable name
g.set_titles("{col_name}")
g.fig.suptitle("Epoch-wise totals")
plt.show()
plt.savefig(os.path.join(figfolder, f"muacontrol_{name}_epochwise.png"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(figfolder, f"muacontrol_{name}_epochwise.pdf"), dpi=300, bbox_inches="tight")

# Create the catplot
g = sns.catplot(data=df_filtered, x='epoch', y='stvalue', col='variable',
                kind='point', col_wrap=4, color='black', height=3, aspect=1.5,
                sharey=False, col_order=variables_to_inspect)
# Adjust the y-axis labels to show the variable name
g.set_titles("{col_name} (standardized)")
g.fig.suptitle("Epoch-wise totals (standardized)")
g.fig.subplots_adjust(top=0.95)
plt.show()
plt.savefig(os.path.join(figfolder, f"muacontrol_{name}_epochwise_standardized.png"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(figfolder, f"muacontrol_{name}_epochwise_standardized.pdf"), dpi=300, bbox_inches="tight")


# Create the catplot
g = sns.catplot(data=df_filtered, x='epoch', y='mmvalue', col='variable',
                kind='point', col_wrap=4, color='black', height=3, aspect=1.5,
                sharey=False, col_order=variables_to_inspect)
# Adjust the y-axis labels to show the variable name
g.set_titles("{col_name} (standardized)")
g.fig.suptitle("Epoch-wise totals (standardized)")
g.fig.subplots_adjust(top=0.95)
plt.show()
plt.savefig(os.path.join(figfolder, f"muacontrol_{name}_epochwise_standardized.png"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(figfolder, f"muacontrol_{name}_epochwise_standardized.pdf"), dpi=300, bbox_inches="tight")

# -----------------------------
# Animal wise version of last three plots (hue)
# -----------------------------

# Create the catplot
g = sns.catplot(data=df_filtered, x='epoch', y='value', col='variable',
                kind='point', col_wrap=4, hue='animal', height=3, aspect=1.5,
                sharey=False, col_order=variables_to_inspect)
# Adjust the y-axis labels to show the variable name
g.set_titles("{col_name}")
g.fig.suptitle("Epoch-wise split by animal")
plt.show()
plt.savefig(os.path.join(figfolder, f"muacontrol_{name}_epochwise_animal.png"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(figfolder, f"muacontrol_{name}_epochwise_animal.pdf"), dpi=300, bbox_inches="tight")

# Create the catplot
g = sns.catplot(data=df_filtered, x='epoch', y='mmvalue', col='variable',
                kind='point', col_wrap=4, height=3, aspect=1.5,
                sharey=False, col_order=variables_to_inspect, hue='animal')
# Adjust the y-axis labels to show the variable name
g.set_titles("{col_name} (minmax)")
g.fig.suptitle("Epoch-wise totals (minmax)")
g.fig.subplots_adjust(top=0.95)
plt.show()
plt.savefig(os.path.join(figfolder, f"muacontrol_{name}_epochwise_minmax_animal.png"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(figfolder, f"muacontrol_{name}_epochwise_standardized_animal.pdf"), dpi=300, bbox_inches="tight")

# Create the catplot
g = sns.catplot(data=df_filtered, x='epoch', y='stvalue', col='variable',
                kind='point', col_wrap=4, height=3, aspect=1.5,
                sharey=False, col_order=variables_to_inspect, hue='animal')
# Adjust the y-axis labels to show the variable name
g.set_titles("{col_name} (standardized)")
g.fig.suptitle("Epoch-wise totals (standardized)")
g.fig.subplots_adjust(top=0.95)
plt.show()
plt.savefig(os.path.join(figfolder, f"{name}_epochwise_standardized_animal.png"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(figfolder, f"{name}_epochwise_standardized_animal.pdf"), dpi=300, bbox_inches="tight")

plt.close('all')

sys.exit(0)
