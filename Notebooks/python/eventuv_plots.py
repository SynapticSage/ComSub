import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler
import os

figfolder = "/Volumes/MATLAB-Drive/Shared/figures/eventuv_python/"
if not os.path.exists(figfolder):
    os.makedirs(figfolder)

# Load the provided CSV file
df = pd.read_parquet('/Volumes/MATLAB-Drive/Shared/figures/tables/eventuv.parquet')
# Display the first few rows of the dataframe
df.head()
# Function to calculate the distance from a point to the line y = x
def distance_to_line(u, v):
    return np.abs(u - v) / np.sqrt(2)
# Function to calculate the Euclidean distance from a point to the origin
def distance_to_origin(u, v):
    return np.sqrt(u**2 + v**2)
# Function to scale a series to the range [0, 1]
def scale(series):
    scaler = MinMaxScaler()
    return scaler.fit_transform(series.values.reshape(-1, 1)).ravel()

def label_epochs(df, time_col='event_time', time_threshold=60*5):
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
df = df.sort_values(['event_time','animal','pattern_cca1','pattern_cca2']).groupby('animal').apply(label_epochs).reset_index()

# Remove rows with NaN in 'event_u_values' or 'event_v_values'
df_clean = df.dropna(subset=['event_u_values', 'event_v_values'])
# Calculate the distance from each point to the line y = x
df_clean['distance_to_line'] = distance_to_line(df_clean['event_u_values'], df_clean['event_v_values'])
# Calculate the Euclidean distance from each point to the origin
df_clean['distance_to_origin'] = distance_to_origin(df_clean['event_u_values'], df_clean['event_v_values'])
# Group the DataFrame and apply the scaling function to the 'distance_to_line' within each group
df_clean['scaled_distance_to_line']   = df_clean.groupby(['genH', 'patterns', 'uv_components', 'animal'])['distance_to_line'].transform(scale)
df_clean['scaled_distance_to_origin'] = df_clean.groupby(['genH', 'patterns', 'uv_components', 'animal'])['distance_to_origin'].transform(scale)
# Calculate 'on_commsub' as one minus the scaled distance to line
df_clean['scaled_on_commsub'] = 1 - df_clean['scaled_distance_to_line']
# Multiply the two distances to get the desired measure
df_clean['scaled_on_commsub_mag'] = df_clean['scaled_on_commsub'] * df_clean['scaled_distance_to_origin']
df_clean['on_commsub'] = 1 - df_clean['distance_to_line']
df_clean['on_commsub_mag'] = df_clean['on_commsub'] * df_clean['distance_to_origin']
df_summary = df_clean.groupby(['genH', 'patterns', 'uv_components', 'animal', 'events']).mean().reset_index()

# Create 'highlow' column based on 'patterns'
df_clean['highlow'] = df_clean['patterns'].apply(lambda x: 'high' if x in [1, 2, 3] else 'low')
# Create a composite column combining 'genH' and 'highlow'
df_clean['genH_highlow'] = df_clean['genH'].astype(str) + "_" + df_clean['highlow']
df_clean['highlow_genH'] = df_clean['highlow'] + "_" + df_clean['genH'].astype(str)
df_clean['magnitude_u'] = np.abs(df_clean['event_u_values'])
df_clean['magnitude_v'] = np.abs(df_clean['event_v_values'])
df_clean[['magnitude_u', 'magnitude_v']].head()

# FILTER
df_clean = df_clean.query('pattern_cca1 == 2 & pattern_cca2 == 7')

# Display the first few rows of the dataframe
df_clean.head()
uv_components = df_clean['uv_components'].unique()
patterns = df_clean['patterns'].unique()
genH_values = df_clean['genH'].unique()

print(df_clean.uv_components.unique())

# PLOT: Hist plot of on/off commsub -------------
# Create a subplot grid with one row for each uv_component and one column for each pattern
def plot_hist(df_clean, thing='on_commsub_mag', pattern_cca2=7):
    fig, axs = plt.subplots(len(uv_components), len(patterns), figsize=(4*len(patterns), 4*len(uv_components)), sharey=True)
    if len(uv_components) == 1:
        axs = [axs]
    if len(patterns) == 1:
        axs = [[ax] for ax in axs]
    q = df_clean[thing].quantile(0.99)
    for i, uv_component in enumerate(uv_components):
        for j, pattern in enumerate(patterns):
            ax = axs[i][j]
            for k, genH in enumerate(genH_values):
                df_sub = df_clean[(df_clean['uv_components'] == uv_component) & 
                                  (df_clean['patterns'] == pattern) & 
                                  (df_clean['genH'] == genH) & 
                                  (df_clean['pattern_cca2'] == pattern_cca2)]
                q_half = df_sub[thing].quantile(0.75)
                # q_90 = df_sub[thing].quantile(0.90)
                q_90 = df_sub[thing].mean()
                g = sns.histplot(df_sub[thing], ax=ax, kde=False, label=genH, binrange=(0, 10))
                palette = sns.color_palette()
                color = palette[k % len(palette)]
                ax.axvline(q_half, color=color, linestyle='--', label=f'75th percentile {genH}')
                ax.axvline(q_90, color=color, linestyle=':')
            if i == 0:
                ax.set_title(f'Pattern {pattern}')
            if j == 0:
                ax.set_ylabel(f'UV Component {int(uv_component)}')
            ax.legend(title='genH')
            ax.set_xlim(0, q)
    plt.tight_layout()
    plt.suptitle(f'{thing}')
    plt.subplots_adjust(top=0.95)
    plt.show()

plt.close('all')
plot_hist(df_clean, thing='on_commsub_mag')
plt.savefig(os.path.join(figfolder, 'plot_hist_thing=on_commsub_mag.png'))
plt.savefig(os.path.join(figfolder, 'plot_hist_thing=on_commsub_mag.pdf'))
plot_hist(df_clean, thing='scaled_on_commsub_mag')
plt.savefig(os.path.join(figfolder, 'plot_hist_thing=scaled_on_commsub_mag.png'))
plt.savefig(os.path.join(figfolder, 'plot_hist_thing=scaled_on_commsub_mag.pdf'))
# plot_hist(df_clean, thing='on_commsub')

# PLOT: Scatter plot of event_u_values vs event_v_values -------------

def plot_hist_aggregated(df, thing='on_commsub_mag', pattern_cca2=7):
    """
    Create a grid of histograms for each combination of pattern, colored by 'genH'.
    Each histogram represents the distribution of the 'thing' column values.
    A vertical dashed line is drawn at the 75th percentile of the 'thing' column for each 'genH'.
    """
    # Identify the unique patterns and genH values
    patterns = df['patterns'].unique()
    genH_values = df['genH'].unique()

    # Create a subplot grid with one column for each pattern
    fig, axs = plt.subplots(1, len(patterns), figsize=(5 * len(patterns), 5), sharey=True)
    
    # Adjust for case when there is only one pattern
    if len(patterns) == 1:
        axs = [axs]

    # Calculate the 99th percentile of the 'thing' column for the x-axis limit
    q = df[thing].quantile(0.99)

    # For each pattern
    for j, pattern in enumerate(patterns):
        ax = axs[j]
        # For each genH value
        for k, genH in enumerate(genH_values):
            # Filter the DataFrame for the current pattern, genH, and pattern_cca2
            df_sub = df[(df['patterns'] == pattern) & (df['genH'] == genH) & (df['pattern_cca2'] == pattern_cca2)]
            # Calculate the 75th percentile of the 'thing' column
            q_half = df_sub[thing].quantile(0.75)
            # Plot a histogram of the 'thing' column values
            g = sns.histplot(df_sub[thing], ax=ax, kde=False, label=genH, binrange=(0, 10), stat='density', common_norm=False)
            # Draw a vertical dashed line at the 75th percentile
            palette = sns.color_palette()
            color = palette[k % len(palette)]
            ax.axvline(q_half, color=color, linestyle='--', label=f'75th percentile {genH}')
        # Set titles and labels
        ax.set_title(f'Pattern {pattern}')
        if j == 0:
            ax.set_ylabel(f'UV Component Aggregated')
        ax.legend(title='genH')
        ax.set_xlim(0, q)

    plt.tight_layout()
    plt.suptitle(f'{thing}')
    plt.subplots_adjust(top=0.95)
    plt.show()

# Test the function with the current DataFrame
plot_hist_aggregated(df_clean, thing='on_commsub_mag')
plt.savefig(os.path.join(figfolder, 'plot_hist_aggregated_thing=on_commsub_mag.png'))
plot_hist_aggregated(df_clean, thing='scaled_on_commsub_mag')
plt.savefig(os.path.join(figfolder, 'plot_hist_aggregated_thing=scaled_on_commsub_mag.png'))


# PLOT: Bar plot of number of unique animals per uv_component, animal, and genH 

# Group by 'uv_component', 'animal', and 'genH' and count unique animals
def plot_animal_counts(df):
    grouped = df.groupby(['uv_components', 'genH']).animal.nunique().reset_index()
    grouped.rename(columns={'animal': 'animal_count'}, inplace=True)

    # Create the bar plot
    plt.figure(figsize=(15, 6))
    sns.barplot(x='uv_components', y='animal_count', hue='genH', data=grouped)

    plt.title('Number of Unique Animals per uv_component, animal, and genH')
    plt.ylabel('Number of Unique Animals')
    plt.xlabel('UV Component')
    plt.legend(title='genH')
    plt.show()

plot_animal_counts(df_clean)

# PLOT: 
import seaborn as sns
import matplotlib.pyplot as plt
# Create a FacetGrid that arranges the data by uv_components (rows) and patterns (columns)
g = sns.FacetGrid(df_clean, col='patterns', row='uv_components', hue='genH',
                  height=4, aspect=1, palette='viridis')
# Apply a scatterplot onto each facet with increased alpha transparency
g.map(sns.scatterplot, 'event_u_values', 'event_v_values', alpha=0.05)
# Add a diagonal line, text, and vertical/horizontal lines at x=0 and y=0 for each subplot
def add_lines(x, y, **kwargs):
    ax = plt.gca()
    ax.axline((0, 0), slope=1, color='black', linestyle='--')
    ax.text(0.6, 0.7, 'comm sub', transform=ax.transAxes, color='black')
    ax.axvline(0, color='black', linestyle=':', linewidth=0.5)
    ax.axvline(5, color='black', linestyle=':', linewidth=0.5)
    ax.axvline(-5, color='black', linestyle=':')
    ax.axhline(-5, color='black', linestyle=':', linewidth=0.5)
    ax.axhline(5, color='black', linestyle=':', linewidth=0.5)
    ax.axhline(0, color='black', linestyle=':')

g.map(add_lines, 'event_u_values', 'event_v_values')
# Add a legend
g.add_legend()
# Display the plot
plt.show()

# ---------------------------------------------
# Calculate the projection score onto the x=y line
df_clean['projection_score'] = (df_clean['event_u_values'] + df_clean['event_v_values']) / np.sqrt(2)
# Calculate the score for distance perpendicular to the x=y line
df_clean['perpendicular_score'] = np.abs(df_clean['event_u_values'] - df_clean['event_v_values']) / np.sqrt(2)
df_clean[['projection_score', 'perpendicular_score']].head()
df_clean['abs_projection_score'] = df_clean['projection_score'].abs()
df_clean['abs_perpendicular_score'] = df_clean['perpendicular_score'].abs()
df_clean['proj_over_perp'] = df_clean['abs_projection_score'] / df_clean['abs_perpendicular_score']
df_clean['abs_proj_over_perp'] = df_clean['abs_projection_score'] / df_clean['abs_perpendicular_score']

def plot_bar(data, y_value, query=None):
    plt.figure(figsize=(12, 6))
    if query:
        data = data.query(query)
    g=sns.barplot(data=data, x='patterns', y=y_value, hue='genH_highlow', errorbar=('ci', 99))
    plt.title(f'Bar Plot for {y_value}\n(99% CI)')
    plt.show()

# Function that takes the sqrt(sum of uv_components squared) 
def mag2(X):
    X = np.sqrt(np.sum(X**2, axis=0))
    return X

# Plotting for projection score
plt.close('all')
index = ['animal','genH', 'genH_highlow', 'highlow', 'patterns', 'event_time']
df_clean_uvgroup = df_clean.groupby(index).agg({'projection_score': mag2,
                                                'perpendicular_score': mag2,
                                                'abs_projection_score': mag2,
                                                'abs_perpendicular_score': mag2,
                                                'proj_over_perp': 'mean',
                                                'abs_proj_over_perp': 'mean'}).reset_index()
plot_bar(df_clean_uvgroup, 'projection_score')
plot_bar(df_clean_uvgroup, 'perpendicular_score')
plot_bar(df_clean, 'abs_projection_score')
plot_bar(df_clean, 'abs_perpendicular_score')
plot_bar(df_clean, 'proj_over_perp')
plot_bar(df_clean, 'abs_proj_over_perp')

def plot_bar_comp(data, y_value):
    sns.catplot(data=data, x='patterns', y=y_value, hue='genH_highlow', errorbar=('ci', 99),
                row='uv_components', col='genH', kind='bar', sharey=True,
                height=2, aspect=1.5)
    plt.suptitle(f'Bar Plot for {y_value}\n(99% CI)')
    plt.show()

# Plotting for projection score
plt.close('all')
plot_bar_comp(df_clean, 'projection_score')
plot_bar_comp(df_clean, 'perpendicular_score')
plot_bar_comp(df_clean, 'abs_projection_score')
plot_bar_comp(df_clean, 'abs_perpendicular_score')

def plot_bar_animal(data, y_value):
    sns.catplot(data=data, x='patterns', y=y_value, hue='highlow_genH',
                errorbar=('ci', 99), row='uv_components', col='animal',
                kind='bar', sharey=False,
                height=1.5, aspect=1.5)
    plt.suptitle(f'Bar Plot for {y_value}\n(99% CI)')
    plt.show()

plt.close('all')
plot_bar_animal(df_clean, 'abs_projection_score')
plot_bar_animal(df_clean, 'abs_perpendicular_score')
plot_bar_animal(df_clean, 'projection_score')
plot_bar_animal(df_clean, 'perpendicular_score')

def plot_bar_anim_summary(data, y_value, comp_max=None):
    if comp_max:
        data = data.query(f'uv_components <= {comp_max}')
    sns.catplot(data=data, x='patterns', y=y_value, hue='genH_highlow',
                errorbar=('ci', 99), col='animal',
                kind='bar', sharey=False,
                height=1.5, aspect=1.5)
    plt.suptitle(f'Bar Plot for {y_value}\n(99% CI)' + \
                 f' for {comp_max} components' if comp_max else '')
    plt.show()

plt.close('all')
plot_bar_anim_summary(df_clean, 'abs_projection_score',    comp_max=5)
plot_bar_anim_summary(df_clean, 'abs_perpendicular_score', comp_max=5)
plot_bar_anim_summary(df_clean, 'projection_score',        comp_max=5)
plot_bar_anim_summary(df_clean, 'perpendicular_score',     comp_max=5)

# ------------
# DIMENSIONALITY REDUCTION
# ------------
index = ['events', 'animal','genH', 'genH_highlow', 'highlow', 'patterns', 'event_time', 'epoch']
# Pivoting the data to get u_values for each uv_component
df_u = df_clean.pivot_table(index=index,
                            columns='uv_components', 
                            values='event_u_values', 
                            aggfunc='first').reset_index()
# Renaming columns for clarity
df_u.columns = [str(col) + '_u' if isinstance(col, float) else col for col in df_u.columns]
# Pivoting the data to get v_values for each uv_component
df_v = df_clean.pivot_table(index=index,
                            columns='uv_components', 
                            values='event_v_values', 
                            aggfunc='first').reset_index()
# Renaming columns for clarity
df_v.columns = [str(col) + '_v' if isinstance(col, float) else col for col in df_v.columns]
# Merging the two dataframes to get the final matrix
df_matrix = df_u.merge(df_v, on=index)
# Removing 'events' and 'animal' columns to get the matrix
matrix = df_matrix.drop(columns=['events', 'animal'])
matrix.head()

# ------------ CCA ------------

from sklearn.cross_decomposition import CCA
from sklearn.preprocessing import StandardScaler

# Choose the first unique animal as the reference
reference_animal = "JS21"
df_matrix = df_matrix.reset_index()

# Split the data into u and v components for CCA
X_ref = df_matrix[df_matrix['animal'] == reference_animal].filter(like='_u')
Y_ref = df_matrix[df_matrix['animal'] == reference_animal].filter(like='_v')

# Fit the CCA model using the reference animal's data
cca = CCA(n_components=X_ref.shape[1])
cca.fit(X_ref, Y_ref)

# Project the data of all animals onto the canonical space of the reference animal
df_matrix_projected = df_matrix.copy()
u_cols = [col for col in df_matrix.columns if '_u' in col]
v_cols = [col for col in df_matrix.columns if '_v' in col]

for animal in df_matrix['animal'].unique():
    print(f'Projecting {animal} onto {reference_animal}')
    X_animal = df_matrix[df_matrix['animal'] == animal][u_cols]
    Y_animal = df_matrix[df_matrix['animal'] == animal][v_cols]
    X_c, Y_c = cca.transform(X_animal, Y_animal)
    df_matrix_projected.loc[df_matrix_projected['animal'] == animal, u_cols] = X_c
    df_matrix_projected.loc[df_matrix_projected['animal'] == animal, v_cols] = Y_c
df_matrix_projected.set_index(index, inplace=True)


# Mean-center each matrix by animal
for animal in df_matrix['animal'].unique():
    print(f'Mean-centering {animal}')
    animal_data = df_matrix_projected.loc[df_matrix_projected.index.get_level_values('animal') == animal]
    mean_centered_data = StandardScaler(with_std=False).fit_transform(animal_data.values)
    df_matrix_projected.loc[df_matrix_projected.index.get_level_values('animal') == animal] = mean_centered_data
df_matrix.set_index(index, inplace=True)

def perform_umap(df_matrix:pd.DataFrame, n_neighbors=15, min_dist=0.1,
                 n_components=2, metric='euclidean'):
    """
    Perform UMAP on the input matrix.
    Parameters:
    - matrix (DataFrame or numpy array): The input data matrix.
    - n_neighbors (int): UMAP parameter. Size of local neighborhood.
    - min_dist (float): UMAP parameter. Minimum distance between points in the low-dimensional representation.
    - n_components (int): Number of dimensions in the low-dimensional space.
    - metric (str): Metric used to compute distances.
    Returns:
    - DataFrame: UMAP-transformed data.
    """
    import umap
    if any([x in df_matrix.columns for x in index]):
        df_matrix = df_matrix.set_index(index)
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, n_components=n_components, metric=metric)
    print(f"Performing UMAP with {n_neighbors} neighbors, {min_dist} min_dist, {n_components} components, and {metric} metric...")
    embedding = reducer.fit_transform(df_matrix.values)
    print(f"...done. UMAP embedding shape: {embedding.shape}")
    return pd.DataFrame(embedding, columns=[f"dim_{i}" for i in range(1, n_components + 1)],
                        index=df_matrix.index)

um =  perform_umap(df_matrix)
um = um.reset_index()
um_anim = perform_umap(df_matrix_projected)
um_anim = um_anim.reset_index()

def plot_umap(matrix:pd.DataFrame, row=None, col=None, hue=None, sample=None,
              col_wrap=None, **kws):
    """
    Plot the UMAP-transformed data.
    Parameters:
    - matrix (DataFrame): UMAP-transformed data.
    - row (str, optional): Column in matrix to use for row-wise subplotting.
    - col (str, optional): Column in matrix to use for column-wise subplotting.
    - hue (str, optional): Column in matrix to use for color encoding.
    Returns:
    - FacetGrid: Seaborn FacetGrid object with scatter plots.
    """
    if sample:
        matrix = matrix.sample(sample)
    g = sns.FacetGrid(matrix, row=row, col=col, hue=hue, height=5, aspect=1,
                      col_wrap=col_wrap)
    g.map(plt.scatter, "dim_1", "dim_2", **kws)
    g.add_legend()
    return g

plot_umap(um, hue='animal', sample=10_000, alpha=0.1)
plot_umap(um_anim, hue='animal', sample=10_000, alpha=0.1)

plt.close('all')
plot_umap(um, hue='genH', col="animal", sample=10_000, alpha=0.1, col_wrap=3)
plot_umap(um_anim, hue='genH', col="animal", sample=10_000, alpha=0.1, col_wrap=3)

plt.close('all')
plot_umap(um, hue='genH_highlow', col="animal", sample=10_000, alpha=0.1, col_wrap=3)
plot_umap(um_anim, hue='genH_highlow', col="animal", sample=10_000, alpha=0.1, col_wrap=3)

plt.close('all')
plot_umap(um, hue='genH_highlow', row="patterns", col="animal", sample=10_000, alpha=0.1)

um =  perform_umap(df_matrix, n_components=3)
um = um.reset_index()
um_anim = perform_umap(df_matrix_projected, n_components=3)
um_anim = um_anim.reset_index()

def plot_umap_3d(matrix, row=None, col=None, hue=None, sample=None,
                 col_wrap=None, **kws):
    """
    Plot the UMAP-transformed data in 3D.
    Parameters:
    - matrix (DataFrame): UMAP-transformed data.
    - row (str, optional): Column in matrix to use for row-wise subplotting.
    - col (str, optional): Column in matrix to use for column-wise subplotting.
    - hue (str, optional): Column in matrix to use for color encoding.
    Returns:
    - FacetGrid: Seaborn FacetGrid object with scatter plots.
    """
    from mpl_toolkits.mplot3d import Axes3D
    if 'magnitude_u' in matrix.columns:
        matrix = matrix.drop(columns=['magnitude_u'])
    if 'magnitude_v' in matrix.columns:
        matrix = matrix.drop(columns=['magnitude_v'])
    if sample:
        matrix = matrix.sample(sample)
    # Create a grid of plots
    g = sns.FacetGrid(matrix, row=row, col=col, hue=hue, height=5, aspect=1,
                      col_wrap=col_wrap, subplot_kws={'projection': '3d'})
    hue_mapping = {unique_value: index for index, unique_value in
                   enumerate(matrix[hue].unique())}
    hue_order = list(hue_mapping.keys())
    for ax, data in g.facet_data():
        palette = sns.utils.get_color_cycle()
        c = g._get_palette(data, hue, hue_order, palette)[ax[-1]]
        print("ax", ax)
        if g.axes.ndim == 1:
            AX = g.axes[ax[1]]
        else:
            AX = g.axes[ax[0], ax[1]]
        AX.scatter(data["dim_1"], data["dim_2"], data["dim_3"], c=c, alpha=0.1)
    # # For each axis object, set it to 3D and plot the 3D scatter plot
    # for ax in g.axes.flat:
    #     ax.remove()
    #     ax = g.fig.add_subplot(ax, projection='3d')
    #     ax.scatter(matrix["dim_1"], matrix["dim_2"], matrix["dim_3"], c=matrix[hue])
    g.add_legend()
    return g

plt.close('all')
plot_umap_3d(um, hue='animal', sample=10_000, alpha=0.1)
plt.savefig(os.path.join(figfolder,'umap_3d_hue=animal.png'), dpi=300)
plt.savefig(os.path.join(figfolder,'umap_3d_hue=animal.pdf'), dpi=300)

plt.close('all')
plot_umap_3d(um, hue='genH', col="animal", sample=10_000, alpha=0.4,
             col_wrap=3)
plt.savefig(os.path.join(figfolder,'umap_3d_hue=genH_col=animal.png'), dpi=300)
plt.savefig(os.path.join(figfolder,'umap_3d_hue=genH_col=animal.pdf'), dpi=300)

plt.close('all')
plot_umap_3d(um, hue='genH_highlow', col="animal", sample=5_000, alpha=1,
             col_wrap=3)
plt.savefig(os.path.join(figfolder,'umap_3d_hue=genH_highlow_col=animal.png'), dpi=300)
plt.savefig(os.path.join(figfolder,'umap_3d_hue=genH_highlow_col=animal.pdf'), dpi=300)

plt.close('all')
plot_umap_3d(um, hue='genH_highlow', row="patterns", col="animal",
             sample=10_000, alpha=0.1)
plt.savefig(os.path.join(figfolder,'umap_3d_hue=genH_highlow_row=patterns_col=animal.png'), dpi=300)
plt.savefig(os.path.join(figfolder,'umap_3d_hue=genH_highlow_row=patterns_col=animal.pdf'), dpi=300)

plot_umap_3d(um_anim, hue='genH_highlow', col="animal", sample=10_000,
alpha=0.9, col_wrap=3)
plt.savefig(os.path.join(figfolder,'umap_3d_rotanim_hue=genH_highlow_col=animal.png'), dpi=300)
plt.savefig(os.path.join(figfolder,'umap_3d_rotanim_hue=genH_highlow_col=animal.pdf'), dpi=300)

# plot_umap_3d(um_anim, hue='animal', sample=10_000, alpha=0.1)
# plot_umap_3d(um_anim, hue='genH', col="animal", sample=10_000, alpha=0.1,
#              col_wrap=3)

# ========
# UV magnitude OVER TIME 󰥔
# ========

# TODO: CUT OUT DIFF GREATER THAN 1 minute

# Calculate magnitude for each u and v component
df_matrix['magnitude_u'] = np.sqrt(df_matrix['1.0_u']**2 +
                                   df_matrix['2.0_u']**2 +
                                   df_matrix['3.0_u']**2 +
                                   df_matrix['4.0_u']**2 +
                                   df_matrix['5.0_u']**2)
df_matrix['magnitude_v'] = np.sqrt(df_matrix['1.0_v']**2 +
                                   df_matrix['2.0_v']**2 +
                                   df_matrix['3.0_v']**2 +
                                   df_matrix['4.0_v']**2 +
                                   df_matrix['5.0_v']**2)

# Reset the index for the operation
df_matrix_reset = df_matrix.reset_index()
# Calculate the linspace for each group
def assign_time(group, time='events'):
    group_size = len(group)
    first_event = group[time].iloc[0]
    last_event = group[time].iloc[-1]
    group['time'] = np.linspace(first_event, last_event, group_size)
    return group
use = 'event_time'
df_matrix_time = df_matrix_reset.groupby(['animal', 'genH']).apply(lambda x:
                                                                   assign_time(x, time=use))
# Setting the new dataframe index back to its original structure
df_matrix_time = df_matrix_time.set_index([use, 'animal', 'genH', 'genH_highlow', 'highlow', 'patterns'])
df_matrix_time.head()


# Subtract minimum time per animal
df_matrix = df_matrix_time.reset_index()
df_matrix['time'] = df_matrix['time'] - df_matrix.groupby('animal')['time'].transform('min')

#9 Melt the data for easy plotting
df_melted = df_matrix.melt(id_vars=['epoch','time'], value_vars=['magnitude_u', 'magnitude_v'], var_name='component', value_name='magnitude')
# Plot using relplot
print("Plotting...")
g = sns.relplot(data=df_melted, x='time', y='magnitude', kind='line', hue='component', height=5, aspect=2, col='epoch', col_wrap=4)
g.set_axis_labels("Event Time", "Magnitude")
g.tight_layout()
plt.show()
plt.close()

# Define window size
window_size = 20
# Apply rolling mean
df_matrix['magnitude_u_smooth'] = df_matrix.groupby(['animal', 'genH'])['magnitude_u'].transform(lambda x: x.rolling(window_size).mean())
df_matrix['magnitude_v_smooth'] = df_matrix.groupby(['animal', 'genH'])['magnitude_v'].transform(lambda x: x.rolling(window_size).mean())

df_melted = df_matrix.melt(id_vars=['time','epoch','animal','genH'], value_vars=['magnitude_u_smooth', 'magnitude_v_smooth'], var_name='component', value_name='magnitude')
# Plot using relplot
print("Plotting...")
g = sns.relplot(data=df_melted, x='time', y='magnitude', kind='line',
                hue='component', height=5, aspect=2, alpha=0.5, col='epoch',
                row='animal', facet_kws=dict(sharex=False))
g.set_axis_labels("Event Time", "Magnitude")
g.tight_layout()
plt.show()
plt.savefig(os.path.join(figfolder,'match_to_genH_overtime_hue=component_row=animal_col=epoch.png'), dpi=300)
plt.savefig(os.path.join(figfolder,'match_to_genH_overtime_hue=component_row=animal_col=epoch.pdf'), dpi=300)

g = sns.relplot(data=df_melted, x='time', y='magnitude', kind='line',
                hue='genH', height=5, aspect=2, alpha=0.5, col='epoch',
                row='animal', facet_kws=dict(sharex=False))
g.set_axis_labels("Event Time", "Magnitude")
g.tight_layout()
plt.show()
plt.savefig(os.path.join(figfolder,'match_to_genH_overtime_hue=component_row=animal_col=epoch.png'), dpi=300)
plt.savefig(os.path.join(figfolder,'match_to_genH_overtime_hue=component_row=animal_col=epoch.pdf'), dpi=300)


# Define number of bins
num_bins = 100
# Bin the event_time data
df_matrix['time_bin'] = df_matrix.groupby(['animal', 'genH', 'epoch'])['time'].transform(lambda x: pd.cut(x, num_bins, labels=range(num_bins)))


# Plot using relplot
print("Plotting...")
g = sns.relplot(data=df_matrix, x='time_bin', y='magnitude_u', kind='line', hue='genH', height=5, aspect=2)
g.set_axis_labels("Time Bin", "Magnitude")
g.tight_layout()
plt.show()


# Plot using relplot
print("Plotting...")
g = sns.relplot(data=df_matrix, x='time_bin', y='magnitude_u', kind='line', hue='genH_highlow', height=5, aspect=2)
g.set_axis_labels("Time Bin", "Magnitude")
g.tight_layout()
plt.show()



# Plot using relplot
print("Plotting...")
g = sns.relplot(data=df_matrix, x='time_bin', y='magnitude_v', kind='line', hue='genH', height=5, aspect=2)
g.set_axis_labels("Time Bin", "Magnitude")
g.tight_layout()
plt.show()


# Plot using relplot
print("Plotting...")
g = sns.relplot(data=df_matrix, x='time_bin', y='magnitude_v', kind='line', hue='genH_highlow', height=5, aspect=2)
g.set_axis_labels("Time Bin", "Magnitude")
g.tight_layout()
plt.show()


# Plot using relplot
print("Plotting...")
g = sns.relplot(data=df_matrix.query('genH == "coherence"'), x='time_bin', y='magnitude_u', kind='line', hue='genH_highlow', height=5, aspect=2)
g.set_axis_labels("Time Bin", "Magnitude")
g.tight_layout()
plt.show()

