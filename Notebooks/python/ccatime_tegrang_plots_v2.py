import scipy.stats
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
from tqdm import tqdm
import os
import imageio
import networkx as nx

animal   = ''
spectype = 'power'
intermediate = 'midpattern=true'
folder   = f'/Volumes/MATLAB-Drive/Shared/figures/{intermediate}/tables/'
plotfolder = f'/Volumes/MATLAB-Drive/Shared/figures/{intermediate}/cca_tegrang'
if not os.path.exists(plotfolder):
    os.makedirs(plotfolder)
color_map = {'Cavg': 'blue', 'S1': 'green', 'S2': 'red', 'wpli': 'purple'}

# Load the data
data = pd.read_csv(os.path.join(folder, f'{animal}{spectype}_granger_causality_directionality.csv'))
data.head()

# --- Transfer entropy comparison ---
# Calculate the proportion of times transfer entropy is greater than surrogate transfer entropy
# Add a new column for the comparison between the original and surrogate transfer entropy
data['te_greater'] = data['transfer_entropy'] > data['surrogate_transfer_entropy']
# Replace the 'direction' column with a more readable version
# data['te_direction'] = data['te_direction'].replace({0: 'column1 -> column2', 1: 'column1 <- column2'})
te_proportion = data.groupby(['column1', 'column2'])['te_greater'].mean()
data = data.set_index(['column1', 'column2']).assign(te_proportion=te_proportion).reset_index()
# Create a pivot table for the heatmap
te_proportion_heatmap_data = te_proportion.reset_index().pivot_table(index='column2', columns='column1', values='te_greater')
# Create a lower diagonal greater than upper diagonal version

def average_lower_diagonal(df):
    df_copy = df.copy()
    lower_diagonal = np.tril(df_copy.values)
    # upper_diagonal = -np.triu(df_copy.values, k=1)
    upper_diagonal = np.triu(df_copy.values)
    averaged_values = (lower_diagonal + upper_diagonal.T) / 2
    averaged_values[averaged_values == 0] = np.nan
    df_averaged = pd.DataFrame(averaged_values, index=df_copy.index, columns=df_copy.columns)
    return df_averaged
te_proportion_heatmap_data_averaged = average_lower_diagonal(te_proportion_heatmap_data)
plt.clf(); sns.heatmap(te_proportion_heatmap_data_averaged, cmap="coolwarm", center=0.5)
plt.subplots_adjust(left=0.3, bottom=0.3)
plt.clf(); sns.heatmap(te_proportion_heatmap_data, cmap="coolwarm", center=0.5)
plt.subplots_adjust(left=0.3, bottom=0.3)

def get_significant_interactions(df, max_lag_order, field="magnitude", pfield="pvalue"):
    df_max_lag = df[df['lag'] == max_lag_order]
    grouped_max_lag = df_max_lag.groupby(['column1', 'column2', 'direction'])

    def combine_and_average(group):
        combined_pstat, combined_pvalue = scipy.stats.combine_pvalues(group[pfield], method='fisher')
        avg_magnitude = group[field].mean()
        if group.direction.iloc[0] == "x->y":
            pass
        elif group.direction.iloc[0] == "y->x":
            avg_magnitude = -avg_magnitude
            combined_pstat = -combined_pstat
        else:
            raise ValueError(f"Invalid direction: {group.direction.iloc[0]}")
        return pd.Series([combined_pstat, combined_pvalue, avg_magnitude, combined_pvalue < 0.05], index=['combined_pstat', 'combined_pvalue', 'avg_magnitude', 'is_significant'])
    
    combined_max_lag = grouped_max_lag.apply(combine_and_average)
    significant_max_lag = combined_max_lag[combined_max_lag['is_significant']]
    heatmap_data = significant_max_lag.pivot_table(index='column1', columns='column2', values='avg_magnitude')
    significant_max_lag['signed_log_magnitude'] = np.sign(significant_max_lag['avg_magnitude']) * np.log10(np.abs(significant_max_lag['avg_magnitude']))
    heatmap_data_signed_log = significant_max_lag.pivot_table(index='column1', columns='column2', values='signed_log_magnitude')
    return significant_max_lag, heatmap_data, heatmap_data_signed_log

# =============================================================================
# Apply the average_lower_diagonal function to the heatmap data for each lag
# =============================================================================
def gen_heatmaps(heatmap_data, heatmap_data_signed_log, tit, average_lower=True):
    if average_lower:
        func = average_lower_diagonal
    else:
        func = lambda x : x
    averaged_heatmap_data = \
            {lag: func(df) for lag, df in heatmap_data.items()}
    averaged_heatmap_data_signed_log = \
            {lag: func(df) for lag, df in heatmap_data_signed_log.items()}
    # -----
    plt.figure()
    sns.heatmap(averaged_heatmap_data[1], cmap="coolwarm", center=0)
    plt.show()
    plt.gcf().suptitle(f'{animal} {tit.replace("_","").capitalize()}')
    # set colorbar label
    plt.gcf().axes[-1].set_ylabel('Magnitude')
    plt.savefig(os.path.join(plotfolder, f'{animal}{spectype}_{tit}_{average_lower}.png'))
    plt.savefig(os.path.join(plotfolder, f'{animal}{spectype}_{tit}_{average_lower}.pdf'))
    # Let's visualize the first lag for each lag
    # Adjust the width of the plot by changing the first value in figsize
    # -----
    fig, ax = plt.subplots(len(averaged_heatmap_data), 1, figsize=(8, 8 * len(averaged_heatmap_data)), squeeze=False)
    for i, data in enumerate(averaged_heatmap_data):
        sns.heatmap(averaged_heatmap_data[data], cmap="coolwarm", center=0, ax=ax[i,0], cbar=False)
        ax[i,0].axis('equal')
        # Adjusting y-axis label position and tick label size
        ax[i,0].set_yticklabels(ax[i,0].get_yticklabels(), fontsize=4)
        ax[i,0].set_xticklabels(ax[i,0].get_xticklabels(), fontsize=4)
    fig.tight_layout()
    fig.suptitle(f'{animal} {tit.replace("_","").capitalize()}')
    fig.savefig(os.path.join(plotfolder, f'{animal}{spectype}_{tit}_lag_{average_lower}.png'))
    fig.savefig(os.path.join(plotfolder, f'{animal}{spectype}_{tit}_lag_{average_lower}.pdf'))
    # -----
    plt.figure()
    sns.heatmap(averaged_heatmap_data_signed_log[1], cmap="coolwarm", center=0)
    plt.show()
    plt.gcf().suptitle(f'{animal} {tit.replace("_","").capitalize()}')
    # set colorbar label
    plt.gcf().axes[-1].set_ylabel('Signed Log Magnitude')
    plt.savefig(os.path.join(plotfolder, f'{animal}{spectype}_{tit}_signed_log_{average_lower}.png'))
    plt.savefig(os.path.join(plotfolder, f'{animal}{spectype}_{tit}_signed_log_{average_lower}.pdf'))
    # Let's visualize the first lag for each lag
    # Adjust the width of the plot by changing the first value in figsize
    fig, ax = plt.subplots(len(averaged_heatmap_data), 1, figsize=(8, 8 * len(averaged_heatmap_data)), squeeze=False)
    for i, data in enumerate(averaged_heatmap_data):
        sns.heatmap(averaged_heatmap_data[data], cmap="coolwarm", center=0, ax=ax[i,0], cbar=False)
        ax[i,0].axis('equal')
        # Adjusting y-axis label position and tick label size
        ax[i,0].set_yticklabels(ax[i,0].get_yticklabels(), fontsize=4)
        ax[i,0].set_xticklabels(ax[i,0].get_xticklabels(), fontsize=4)
    fig.tight_layout()
    fig.suptitle(f'{animal} {tit.replace("_","").capitalize()}')
    fig.savefig(os.path.join(plotfolder, f'{animal}{spectype}_{tit}_lag_{average_lower}.png'))
    fig.savefig(os.path.join(plotfolder, f'{animal}{spectype}_{tit}_lag_{average_lower}.pdf'))
    # --
    def create_heatmap(df, lag, vmin, vmax, cmap="coolwarm"):
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.heatmap(df, cmap=cmap, center=0, vmin=vmin, vmax=vmax)
        plt.title(f"Lag order: {lag}")
        plt.savefig(os.path.join(plotfolder,f"heatmap_lag_{lag}_{average_lower}.png"))
        plt.close()
    def get_global_vals(thing):
        global_min = min(df.min().min() for df in thing.values())
        global_max = max(df.max().max() for df in thing.values())
        return global_min, global_max
    global_min, global_max = get_global_vals(averaged_heatmap_data_signed_log)
    # --
    # Generate and save heatmaps
    for lag, df in averaged_heatmap_data_signed_log.items():
        create_heatmap(df, lag, global_min, global_max)
    # --
    # Create GIF
    images = []
    for lag in range(1, lag_max + 1):
        images.append(imageio.imread(os.path.join(plotfolder,f"heatmap_lag_{lag}_{average_lower}.png")))
    # Append the first image to make the GIF loop
    images.append(images[0])
    imageio.mimsave(os.path.join(plotfolder,f'{animal}{tit}_log_heatmap_{average_lower}.gif'),
                    images, duration=1, loop=1000)  # duration is the time between frames in seconds
    # Generate and save heatmaps
    global_min, global_max = get_global_vals(averaged_heatmap_data)
    for lag, df in averaged_heatmap_data.items():
        create_heatmap(df, lag, global_min, global_max)
    # --
    # Create GIF
    images = []
    for lag in range(1, lag_max + 1):
        images.append(imageio.imread(os.path.join(plotfolder,f"heatmap_lag_{lag}_{average_lower}.png")))
    # Append the first image to make the GIF loop
    images.append(images[0])
    imageio.mimsave(os.path.join(plotfolder,f'{animal}{tit}_heatmap_{average_lower}.gif'),
                    images, duration=1, loop=1000)  # duration is the time between frames in seconds

# ================================================================
# GRANGER CAUSALITY DIRECTION
# ================================================================
lag_max = data['lag'].max()
significant_max_lag, heatmap_data, heatmap_data_signed_log = {}, {}, {}
for lag in range(1, lag_max+1):
    significant_max_lag[lag], heatmap_data[lag], heatmap_data_signed_log[lag] \
            = get_significant_interactions(data, max_lag_order=lag)

gen_heatmaps(heatmap_data, heatmap_data_signed_log, tit="Granger Causality")
gen_heatmaps(heatmap_data, heatmap_data_signed_log, tit="Granger Causality", average_lower=False)

plt.close('all')

# ================================================================
# TRANSFER ENTROPY
# ================================================================

lag_max = data['lag'].max()
significant_max_lag, heatmap_data, heatmap_data_signed_log = {}, {}, {}
for lag in range(1, lag_max+1):
    significant_max_lag[lag], heatmap_data[lag], heatmap_data_signed_log[lag] \
            = get_significant_interactions(data, max_lag_order=lag, 
                                           field="transfer_entropy", 
                                           pfield="te_greater")

gen_heatmaps(heatmap_data, heatmap_data_signed_log, tit="transfer entropy")
gen_heatmaps(heatmap_data, heatmap_data_signed_log, tit="transfer entropy", average_lower=False)

plt.close('all')
# --- Transfer entropy comparison with signed log ---
# Create a signed log version of the transfer entropy comparison heatmap data
te_proportion_heatmap_data_signed_log = te_proportion_heatmap_data.applymap(signed_log10)

# Visualize the heatmap
sns.heatmap(te_proportion_heatmap_data_signed_log, cmap="coolwarm", center=0)
plt.show()

# =============================================================================
# the following code should be kept the same as it is performing some pre-processing and doesn't have any computation involving TE or surrogate TE.

# --- Data preparation ---
# # Add a new column for the absolute difference between the original and surrogate transfer entropy
# df['te_difference'] = (df['te'] - df['surrogate_te']).abs()
#
# # Define a new function for signed log transformation
# def signed_log10(x):
#     if x == 0:
#         return 0
#     else:
#         return np.sign(x) * np.log10(np.abs(x))
#
# # Apply the function to the 'te_difference' column
# df['te_difference_signed_log'] = df['te_difference'].apply(signed_log10)
#
# # Replace the 'direction' column with a more readable version
# df['direction'] = df['direction'].replace({0: 'column1 -> column2', 1: 'column1 <- column2'})
#
# # Create a pivot table for the signed log difference between the original and surrogate transfer entropy
# te_difference_heatmap_data = df.pivot_table(index='column1', columns='column2', values='te_difference_signed_log')
#
#
# --- Data preparation ---
# Create a pivot table for the comparison between the original and surrogate transfer entropy
te_comparison_heatmap_data = df.pivot_table(index='column1', columns='column2', values='te_greater')


# =============================================================================

# --- Heatmap ---
# Create the heatmap using seaborn
plt.figure(figsize=(10, 8))
sns.heatmap(te_comparison_heatmap_data, cmap='RdYlGn', center=0, annot=True, fmt=".2f", linewidths=.5)
plt.title('Transfer Entropy Comparison (TE > Surrogate TE)', fontsize=15)
plt.xlabel('column2', fontsize=12)
plt.ylabel('column1', fontsize=12)
plt.show()

# =============================================================================

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def generate_heatmap(df, column1, column2, te_col, surrogate_te_col):
    """
    Function to generate a heatmap summarizing the proportion of values where 
    X->Y Transfer Entropy is greater than the surrogate for that pair.

    Parameters:
    - df: DataFrame containing the data
    - column1, column2: Names of the two columns defining the pairs
    - te_col: Name of the column containing the Transfer Entropy values
    - surrogate_te_col: Name of the column containing the surrogate Transfer Entropy values

    Returns:
    - A heatmap plot
    """

    # Initialize an empty DataFrame to hold the proportions
    proportions_df = pd.DataFrame()

    # Get the unique values in column1 and column2
    unique_values = list(set(df[column1].unique()).union(set(df[column2].unique())))

    # Loop over all unique values
    for i, value1 in enumerate(unique_values):
        for j, value2 in enumerate(unique_values):
            if i < j:
                # For each pair, calculate the proportion of TE values greater than the surrogate
                df_subset = df[(df[column1] == value1) & (df[column2] == value2)]
                proportion = (df_subset[te_col] > df_subset[surrogate_te_col]).mean()
                proportions_df.loc[value1, value2] = proportion

    # Create the heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(proportions_df, cmap='coolwarm', center=0.5, annot=True)

    # Labels and title
    plt.title("Proportions of X->Y TE values greater than surrogate")
    plt.xlabel(column1)
    plt.ylabel(column2)

    return plt

# Example usage (replace 'column1', 'column2', 'te', 'surrogate_te' with your actual column names)
# generate_heatmap(df, 'column1', 'column2', 'te', 'surrogate_te')

# =1============================================================================


def generate_proportion_stripplot(df, column1, column2, te_col, surrogate_te_col):
    # Create a new DataFrame to store the proportions
    proportions_df = pd.DataFrame(columns=[column1, column2, 'Proportion'])

    # Iterate over each unique pair of values in column1 and column2
    for val1 in df[column1].unique():
        for val2 in df[column2].unique():
            # Compute the proportion of TE values greater than the surrogate for this pair
            pair_df = df[(df[column1] == val1) & (df[column2] == val2)]
            proportion = (pair_df[te_col] > pair_df[surrogate_te_col]).mean()
            proportions_df = proportions_df.append({column1: val1, column2: val2, 'Proportion': proportion}, ignore_index=True)

    # Generate a stripplot of the proportions
    plt.figure(figsize=(10, 6))
    sns.stripplot(x=column1, y='Proportion', hue=column2, data=proportions_df, jitter=True)
    plt.title(f'Proportion of {te_col} values greater than {surrogate_te_col}')
    plt.ylabel('Proportion')
    plt.show()

# =12============================================================================


def generate_proportion_boxplot(df, column1, column2, te_col, surrogate_te_col):
    # Create a new DataFrame to store the proportions
    proportions_df = pd.DataFrame(columns=[column1, column2, 'Proportion'])

    # Iterate over each unique pair of values in column1 and column2
    for val1 in df[column1].unique():
        for val2 in df[column2].unique():
            # Compute the proportion of TE values greater than the surrogate for this pair
            pair_df = df[(df[column1] == val1) & (df[column2] == val2)]
            proportion = (pair_df[te_col] > pair_df[surrogate_te_col]).mean()
            proportions_df = proportions_df.append({column1: val1, column2: val2, 'Proportion': proportion}, ignore_index=True)

    # Generate a boxplot of the proportions
    plt.figure(figsize=(10, 6))
    sns.boxplot(x=column1, y='Proportion', hue=column2, data=proportions_df)
    plt.title(f'Proportion of {te_col} values greater than {surrogate_te_col}')
    plt.ylabel('Proportion')
    plt.show()

# =1123============================================================================


import networkx as nx

def generate_proportion_network(df, column1, column2, te_col, surrogate_te_col):
    # Create a new network graph
    G = nx.Graph()

    # Add nodes to the graph
    nodes = set(df[column1].unique()).union(set(df[column2].unique()))
    G.add_nodes_from(nodes)

    # Add edges to the graph
    for val1 in df[column1].unique():
        for val2 in df[column2].unique():
            # Compute the proportion of TE values greater than the surrogate for this pair
            pair_df = df[(df[column1] == val1) & (df[column2] == val2)]
            proportion = (pair_df[te_col] > pair_df[surrogate_te_col]).mean()

            # Add an edge to the graph with this proportion as the weight
            G.add_edge(val1, val2, weight=proportion)

    # Draw the graph
    pos = nx.spring_layout(G)  # positions for all nodes

    # nodes
    nx.draw_networkx_nodes(G, pos, node_size=700)

    # edges
    nx.draw_networkx_edges(G, pos, width=[d['weight'] for (u, v, d) in G.edges(data=True)])

    # labels
    nx.draw_networkx_labels(G, pos, font_size=20, font_family='sans-serif')

    plt.title(f'Network of {te_col} values greater than {surrogate_te_col}')
    plt.axis('off')
    plt.show()


import networkx as nx
import matplotlib.pyplot as plt

def generate_graph_from_df(df, pvalue_field='pvalue', magnitude_field='F', significance_level=0.05):
    """
    Generate a directed graph based on Granger causality results.
    
    Parameters:
    - df: DataFrame containing Granger causality results.
    - pvalue_field: column name for p-values to determine edges.
    - magnitude_field: column name for edge weights.
    - significance_level: significance level for p-value to consider an edge.
    
    Returns:
    - G: A networkx DiGraph.
    """
    # Create a directed graph
    G = nx.DiGraph()
    
    # Iterate through rows of the dataframe
    for _, row in df.iterrows():
        source = row['column1']
        target = row['column2']
        
        # Check if the p-value is less than the significance level
        if row[pvalue_field] < significance_level:
            # Add an edge with weight based on the magnitude_field
            G.add_edge(source, target, weight=row[magnitude_field])
    
    return G


# Visualization
def visualize_graph(G, scale=0.001, zero_id=True):
    pos = nx.spring_layout(G, iterations=20_000)
    weights = [0 if (u==v and zero_id) else G[u][v]['weight']*scale for u, v in G.edges()]
    nx.draw(G, pos, with_labels=True, width=weights, node_size=5,
            node_color="skyblue", font_size=15, font_color="red", 
            edge_cmap=plt.cm.Blues)
    plt.title("Granger Causality Network")
    plt.show()
    return G

# Generate a graph using the function
plt.close('all')
G = generate_graph_from_df(data)  # Assuming results_df is your Granger causality results dataframe
G = visualize_graph(G)
plt.savefig(os.path.join(plotfolder, "granger_causality_network.png"), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(plotfolder, "granger_causality_network.pdf"), dpi=300, bbox_inches='tight')

plt.close('all')
for lag in sorted(data.lag.unique()):  # Assuming data is your time series data
    plt.figure()
    G = generate_graph_from_df(data.query(f'lag == {lag}'))
    G = visualize_graph(G)
    plt.pause(0.5)
    plt.suptitle(f'Granger Causality Network (lag={lag})')
    plt.savefig(os.path.join(plotfolder, f"granger_causality_network_lag_{lag}.png"), dpi=300, bbox_inches='tight')



# ================== GRAPHS WITH STATIC NODES ===================================

def generate_graph_from_df(df, pvalue_field='pvalue', magnitude_field='F', significance_level=0.05):
    """
    Generate a directed graph based on Granger causality results.
    
    Parameters:
    - df: DataFrame containing Granger causality results.
    - pvalue_field: column name for p-values to determine edges.
    - magnitude_field: column name for edge weights.
    - significance_level: significance level for p-value to consider an edge.
    
    Returns:
    - G: A networkx DiGraph.
    """
    G = nx.DiGraph()
    for _, row in df.iterrows():
        source = row['column1']
        target = row['column2']
        if row[pvalue_field] < significance_level:
            G.add_edge(source, target, weight=row[magnitude_field])
    return G

def visualize_graph(G, pos, scale=0.001, zero_id=True):
    weights = [0 if (u == v and zero_id) else G[u][v]['weight'] * scale for u, v in G.edges()]
    nx.draw(G, pos, with_labels=True, width=weights, node_size=5, 
            node_color="skyblue", font_size=15, font_color="red", 
            edge_cmap=plt.cm.Blues)
    plt.title("Granger Causality Network")
    plt.show()

# Compute the overall node positions
full_graph = generate_graph_from_df(data)
positions = nx.spring_layout(full_graph, iterations=20_000)

# Visualize the overall graph
plt.figure()
visualize_graph(full_graph, positions)

# Visualize lag-specific graphs
mng = plt.get_current_fig_manager()
plt.figure()
sz = mng.window.screen().size()
for lag in sorted(data.lag.unique()):
    plt.clf()
    lag_graph = generate_graph_from_df(data.query(f'lag == {lag}'))
    visualize_graph(lag_graph, positions)
    mng.resize(sz.width(), sz.height())
    plt.suptitle(f'Granger Causality Network (lag={lag})')
    plt.savefig(os.path.join(plotfolder, f"granger_causality_network_lag_{lag}.png"), dpi=300, bbox_inches='tight')
    plt.pause(0.5)
