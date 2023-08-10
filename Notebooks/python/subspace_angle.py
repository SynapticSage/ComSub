#!/usr/bin/env python3

import argparse
import os
from matplotlib import pyplot as plt
import hdf5storage
import seaborn as sns
import numpy as np
import pandas as pd
from matplotlib import colors
import networkx as nx
import os

# Argument parser
parser = argparse.ArgumentParser(description="Load and process .mat file for neural data analysis.")
parser.add_argument('--file', type=str, default="/Volumes/MATLAB-Drive/Shared/figures/subspaceAngle/type=wpli/subspaceDist_K=6_norm=1_measure=princ_unitvar.mat", help="Path to the .mat file.")
args = parser.parse_args()

# Load the .mat file
file = os.path.expanduser(args.file)
file = os.path.abspath(file)
data = hdf5storage.loadmat(file)
folder, basename = os.path.split(file)
base, ext  = os.path.splitext(basename)
figfolder = os.path.join(folder, base)
if not os.path.exists(figfolder):
    os.makedirs(figfolder)
plt.ion()

if __name__ == "__main__":

    if 'frobDist' not in data:
        print("Frob dist not in data")
        os._exit(0)
    else:
        print(f"PLOTTING {args.file}")

    # Extract the variables
    subspaceDist = data['subspaceDist']
    frobDist     = data['frobDist']
    rowVar       = data['rowVar']


    # Extract the variables and their shapes from the v6 data
    keys = list(data.keys())
    subspaceDist_shape = data['subspaceDist'].shape if 'subspaceDist' in keys else None
    frobDist_shape = data['frobDist'].shape if 'frobDist' in keys else None
    rowVar_shape = data['rowVar'].shape if 'rowVar' in keys else None
    rowVar_content = data['rowVar'] if 'rowVar' in keys else None

    keys, subspaceDist_shape, frobDist_shape, rowVar_shape, rowVar_content

    # Extract rowVar for tick labels
    tick_labels = [str(x[0]) for x in data['rowVar'][0]]

    mean_subspaceDist   = np.mean(subspaceDist, axis=0).squeeze()
    median_subspaceDist = np.median(subspaceDist, axis=0).squeeze()
    mean_frobDist       = np.mean(frobDist, axis=0).squeeze()
    median_frobDist     = np.median(frobDist, axis=0).squeeze()

    # ----------------------------

    # Define colormaps
    subspace_cmap = 'plasma'
    frob_cmap     = 'inferno'

    def plot_heatmaps(scale='linear'):
        # Define range for subspaceDist colormap
        vmin_subspace, vmax_subspace = 0, np.pi/2

        # Plotting with modifications
        fig, axs = plt.subplots(2, 2, figsize=(14, 12))
        if scale == 'linear':
            mean_frobDist       = globals()['mean_frobDist']
            median_frobDist     = globals()['median_frobDist']
            mean_subspaceDist   = globals()['mean_subspaceDist']
            median_subspaceDist = globals()['median_subspaceDist']
        elif scale == 'log':
            mean_frobDist       = np.log(globals()['mean_frobDist'])
            median_frobDist     = np.log(globals()['median_frobDist'])
            mean_subspaceDist   = np.log(globals()['mean_subspaceDist'])
            median_subspaceDist = np.log(globals()['median_subspaceDist'])
        else:
            raise ValueError("Invalid scale value.")

        # Mean subspaceDist heatmap
        cax1 = axs[0, 0].matshow(mean_subspaceDist, cmap=subspace_cmap, aspect='auto', vmin=vmin_subspace, vmax=vmax_subspace)
        cbar1 = fig.colorbar(cax1, ax=axs[0, 0], ticks=[0, np.pi/4, np.pi/2])
        cbar1.ax.set_yticklabels(['0', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$'])  # LaTeX format
        axs[0, 0].set_title("Mean subspaceDist")
        axs[0, 0].set_xticks(np.arange(len(tick_labels)))
        axs[0, 0].set_yticks(np.arange(len(tick_labels)))
        axs[0, 0].set_xticklabels(tick_labels, rotation=45, ha="right", fontsize=10)
        axs[0, 0].set_yticklabels(tick_labels, fontsize=10)

        # Median subspaceDist heatmap
        cax2 = axs[0, 1].matshow(median_subspaceDist, cmap=subspace_cmap, aspect='auto', vmin=vmin_subspace, vmax=vmax_subspace)
        cbar2 = fig.colorbar(cax2, ax=axs[0, 1], ticks=[0, np.pi/4, np.pi/2])
        cbar2.ax.set_yticklabels(['0', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$'])  # LaTeX format
        axs[0, 1].set_title("Median subspaceDist")
        axs[0, 1].set_xticks(np.arange(len(tick_labels)))
        axs[0, 1].set_yticks(np.arange(len(tick_labels)))
        axs[0, 1].set_xticklabels(tick_labels, rotation=45, ha="right", fontsize=10)
        axs[0, 1].set_yticklabels(tick_labels, fontsize=10)

        # Mean frobDist heatmap
        cax3 = axs[1, 0].matshow(mean_frobDist, cmap=frob_cmap, aspect='auto')
        cbar3 = fig.colorbar(cax3, ax=axs[1, 0])
        axs[1, 0].set_title("Mean frobDist")
        axs[1, 0].set_xticks(np.arange(len(tick_labels)))
        axs[1, 0].set_yticks(np.arange(len(tick_labels)))
        axs[1, 0].set_xticklabels(tick_labels, rotation=45, ha="right", fontsize=10)
        axs[1, 0].set_yticklabels(tick_labels, fontsize=10)

        # Median frobDist heatmap
        cax4 = axs[1, 1].matshow(median_frobDist, cmap=frob_cmap, aspect='auto')
        cbar4 = fig.colorbar(cax4, ax=axs[1, 1])
        axs[1, 1].set_title("Median frobDist")
        axs[1, 1].set_xticks(np.arange(len(tick_labels)))
        axs[1, 1].set_yticks(np.arange(len(tick_labels)))
        axs[1, 1].set_xticklabels(tick_labels, rotation=45, ha="right", fontsize=10)
        axs[1, 1].set_yticklabels(tick_labels, fontsize=10)

        plt.tight_layout()
        plt.show()

    plot_heatmaps(scale='linear')
    plt.savefig(os.path.join(figfolder, 'subspaceDist_frobDist.png'), dpi=300)

    plot_heatmaps(scale='log')
    plt.savefig(os.path.join(figfolder, 'subspaceDist_frobDist_log.png'), dpi=300)

    # ----------------------------


    # Helper function to create dataframe for given distance tensor
    def create_dataframe(distance_tensor, row_vars):
        num_samples, num_conditions, _ = distance_tensor.shape
        isample_list, rowvarX_list, rowvarY_list, distance_list = [], [], [], []
        
        for sample in range(num_samples):
            for i in range(num_conditions):
                for j in range(num_conditions):
                    isample_list.append(sample)
                    rowvarX_list.append(row_vars[i])
                    rowvarY_list.append(row_vars[j])
                    distance_list.append(distance_tensor[sample, i, j])
                    
        df = pd.DataFrame({
            'isample': isample_list,
            'rowvarX': rowvarX_list,
            'rowvarY': rowvarY_list,
            'distance': distance_list
        })
        
        return df

    # ----------------------------

    # Create dataframes
    subspaceDf = create_dataframe(data['subspaceDist'], tick_labels)
    frobDf     = create_dataframe(data['frobDist'], tick_labels)

    # ----------------------------

    from scipy.cluster.hierarchy import linkage, leaves_list

    def compute_dendrogram_ordering(matrix):
        """
        Compute the ordering of rows and columns based on hierarchical clustering.
        
        Args:
        - matrix (pd.DataFrame): The input matrix.
        
        Returns:
        - Ordered row indices and ordered column indices.
        """
        # Compute the linkage for rows and columns
        row_linkage = linkage(matrix, method="average")
        col_linkage = linkage(matrix.T, method="average")
        
        # Extract the ordering of rows and columns based on the linkage
        row_order = leaves_list(row_linkage)
        col_order = leaves_list(col_linkage)
        
        return row_order, col_order

    def create_side_by_side_clustermap(df, set_minoffdiag=False):
        """
        Create two side-by-side Seaborn clustermaps (mean and standard error) for a given distance dataframe.
        
        Args:
        - df (pd.DataFrame): The dataframe containing the distance information.
        
        Returns:
        - Two Seaborn clustermap visualizations side by side.
        """
        # Pivot the dataframe to create matrices for mean and standard error
        mean_matrix = df.pivot_table(index="rowvarX", columns="rowvarY", values="distance", aggfunc="mean")
        stderr_matrix = df.pivot_table(index="rowvarX", columns="rowvarY", values="distance", aggfunc=lambda x: np.std(x) / np.sqrt(len(x)))

        # Compute dendrogram ordering for rows and columns
        row_order, col_order = compute_dendrogram_ordering(mean_matrix)
        
        # Reorder the matrices based on the dendrogram ordering
        mean_matrix   = mean_matrix.iloc[row_order, col_order]
        stderr_matrix = stderr_matrix.iloc[row_order, col_order]
        
        # Create the clustermaps with precomputed ordering
        fig, axs = plt.subplots(1, 2, figsize=(20, 8))

        vmax = np.max(np.max(mean_matrix).values)
        # min of non-diagonal elements
        if set_minoffdiag:
            I = np.eye(mean_matrix.shape[0], dtype=bool)
            vmin = np.min(mean_matrix.values[~I])
            mean_matrix.values[I] = np.nan
        else:
            vmin = 0
        
        sns.heatmap(mean_matrix, cmap="viridis", ax=axs[0], cbar_kws={'label': 'Mean Distance'}, vmax=vmax, vmin=vmin)
        axs[0].set_title("Mean Distance")
        
        sns.heatmap(stderr_matrix, cmap="viridis", ax=axs[1], cbar_kws={'label': 'Standard Error'})
        axs[1].set_title("Standard Error")

        plt.tight_layout()
        plt.show()

    # Create side-by-side clustermaps for the subspaceDf with the corrected colorbar for mean distances
    plt.close('all')
    create_side_by_side_clustermap(subspaceDf, set_minoffdiag=True)
    plt.savefig(os.path.join(figfolder, 'subspaceDist_clustermap.png'), dpi=300)
    plt.savefig(os.path.join(figfolder, 'subspaceDist_clustermap.pdf'), dpi=300)
    create_side_by_side_clustermap(frobDf, set_minoffdiag=True)
    plt.savefig(os.path.join(figfolder, 'frobDist_clustermap.png'), dpi=300)
    plt.savefig(os.path.join(figfolder, 'frobDist_clustermap.pdf'), dpi=300)


    
    # ----------------------------
    def create_ordered_barplot_with_options(df, color_col=None, subtract_mean=False, percent_scale=False):
        """
        Create an ordered Seaborn barplot for a given distance dataframe with options for mean subtraction and percent scaling.
        
        Args:
        - df (pd.DataFrame): The dataframe containing the distance information.
        - color_col (str, optional): Column name to color the bars by. Defaults to None.
        - subtract_mean (bool): If True, subtracts the mean of the measurements after NaNing identity elements in the dataframe.
        - percent_scale (bool): If True, scales distances as a percent of the maximum value.
        
        Returns:
        - Seaborn barplot visualization.
        """
        # Pivot the dataframe to create a matrix for mean values
        mean_matrix = df.pivot_table(index="rowvarX", columns="rowvarY", values="distance", aggfunc="mean")
        
        # Compute dendrogram ordering for rows and columns
        row_order, col_order = compute_dendrogram_ordering(mean_matrix.fillna(0))
        
        # Create a new column combining rowvarX and rowvarY to represent each bar
        df["pairs"] = df["rowvarX"] + " - " + df["rowvarY"]
        
        # Subtract the mean and NaN the diagonal if subtract_mean is True
        if subtract_mean:
            df["distance"] = df["distance"] - df["distance"].mean()
            df.loc[df["rowvarX"] == df["rowvarY"], "distance"] = np.nan
        
        # Scale to percent of max value if percent_scale is True
        if percent_scale:
            max_value = df["distance"].max()
            df["distance"] = (df["distance"] / max_value) * 100
        
        # Order the bars based on the dendrogram ordering
        ordered_pairs = [mean_matrix.index[i] + " - " + mean_matrix.columns[j] for i in row_order for j in col_order]
        
        # Create the barplot
        plt.figure(figsize=(20, 10))
        sns.barplot(data=df, x="pairs", y="distance", order=ordered_pairs, hue=color_col, palette="viridis")
        plt.xticks(rotation=90)
        plt.ylabel("Distance (% of Max)" if percent_scale else "Distance")
        plt.title("Ordered Barplot of Distances with Mean Subtraction" if subtract_mean else "Ordered Barplot of Distances")
        plt.tight_layout()
        plt.show()

    # Create the ordered barplot for the frobDf with mean subtraction and percent scaling
    create_ordered_barplot_with_options(frobDf, subtract_mean=True, percent_scale=True)

    # ----------------------------

    # Plot histograms for subspaceDf and frobDf to visualize the distribution of distances

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    # Histogram for subspaceDf
    sns.histplot(subspaceDf["distance"].dropna(), ax=axes[0], bins=30, kde=True, color="blue")
    axes[0].set_title("Distribution of Subspace Distances")
    axes[0].set_xlabel("Subspace Distance")
    axes[0].set_ylabel("Frequency")

    # Histogram for frobDf
    sns.histplot(frobDf["distance"].dropna(), ax=axes[1], bins=30, kde=True, color="green")
    axes[1].set_title("Distribution of Frobenius Distances")
    axes[1].set_xlabel("Frobenius Distance")
    axes[1].set_ylabel("Frequency")

    plt.tight_layout()
    plt.show()
    plt.savefig(os.path.join(figfolder, 'distance_value_distribution.png'), dpi=300)
    plt.savefig(os.path.join(figfolder, 'distance_value_distribution.pdf'))


    # ----------------------------
    # def pairwise_scatter_corrected(df):
    #     import itertools
    #     unique_conditions = df['rowvarX'].unique()
    #     pairs = list(itertools.combinations(unique_conditions, 2))
    #     
    #     # Check if there are pairs to plot
    #     if len(pairs) == 0:
    #         print("Not enough unique conditions to generate pairwise scatter plots.")
    #         return
    #
    #     fig, axes = plt.subplots(len(pairs), figsize=(15, 5*len(pairs)))
    #
    #     for ax, (cond1, cond2) in zip(axes, pairs):
    #         distances_cond1 = df[df['rowvarX'] == cond1]['distance']
    #         distances_cond2 = df[df['rowvarX'] == cond2]['distance']
    #         
    #         sns.scatterplot(x=distances_cond1, y=distances_cond2, ax=ax)
    #         ax.set_xlabel(cond1)
    #         ax.set_ylabel(cond2)
    #         ax.set_title(f"{cond1} vs {cond2}")
    #
    #     plt.tight_layout()
    #     plt.show()
    #
    # # Running the function for the frobDf dataframe
    # pairwise_scatter_corrected(frobDf)

    # ----------------------------


    def plot_network_graph(df, scale=1.5):
    # Create a new graph from the dataframe
        G = nx.Graph()

        for index, row in df.iterrows():
            if row['rowvarX'] == row['rowvarY'] or row['distance'] == 0:
                weight = 0
            else:
                weight = 1/row['distance']  # inverse weight for visualization
            G.add_edge(row['rowvarX'], row['rowvarY'], weight=weight)  # inverse weight for visualization

        pos = nx.spring_layout(G, k=0.5, iterations=1000)  # k = distance between nodes, iterations = number of iterations to run the algorithm
        weights = nx.get_edge_attributes(G, 'weight').values()
        weights = [scale * w for w in weights]  # scale edge weights for visualization

        print("Edge weights:")
        print(weights)

        # Draw the network graph
        plt.figure(figsize=(10, 10))
        nx.draw_networkx_nodes(G, pos, node_size=1000)
        nx.draw_networkx_labels(G, pos)
        nx.draw_networkx_edges(G, pos, width=list(weights))
        
        # Draw edge labels
        edge_labels = nx.get_edge_attributes(G, 'weight')
        edge_labels = {(u, v): '{:.3f}'.format(weight) for (u, v), weight in nx.get_edge_attributes(G, 'weight').items()}
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, label_pos=0.5, font_size=10, font_color='black')

        plt.title('Network Graph of Distances')
        plt.show()
        return G, weights


    # Running the function for the subspaceDf dataframe
    G = plot_network_graph(subspaceDf)
    plt.savefig(os.path.join(figfolder, "graph_subspace"))
    G = plot_network_graph(frobDf)
    plt.savefig(os.path.join(figfolder, "graph_frobdist"))
