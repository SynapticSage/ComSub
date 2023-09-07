import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler
import os
from statsmodels.api import qqplot_2samples
from tqdm import tqdm
tqdm.pandas()

match_mode = False
zscore = True
intermediate = "midpattern=true"
figfolder = f"/Volumes/MATLAB-Drive/Shared/figures/{intermediate}/eventuv_python/"
origin=f'/Volumes/MATLAB-Drive/Shared/figures/{intermediate}/tables/eventuv.parquet'
if not os.path.exists(figfolder):
    os.makedirs(figfolder)

def prep_uv_melt(df_clean):
    # ------------
    # DIMENSIONALITY REDUCTION
    # ------------
    index = ['events', 'animal','genH', 'genH_highlow', 'highlow', 'patterns',
             'event_time', 'epoch', 'highlow_genH', 'pattern_class', 'pattern_cca1', 'pattern_cca2']
    # Pivoting the data to get u_values for each uv_component
    df_u = df_clean.pivot_table(index=index,
                                columns='uv_components', 
                                values='event_u_values', 
                                aggfunc='mean').reset_index()
    # Renaming columns for clarity
    df_u.columns = [str(col) + '_u' if isinstance(col, float) else col for col in df_u.columns]
    # Pivoting the data to get v_values for each uv_component
    df_v = df_clean.pivot_table(index=index,
                                columns='uv_components', 
                                values='event_v_values', 
                                aggfunc='mean').reset_index()
    # Renaming columns for clarity
    df_v.columns = [str(col) + '_v' if isinstance(col, float) else col for col in df_v.columns]
    # Merging the two dataframes to get the final matrix
    df_matrix = df_u.merge(df_v, on=index)
    # Removing 'events' and 'animal' columns to get the matrix
    matrix = df_matrix.drop(columns=['events', 'animal'])
    matrix.head()
    # drop all rows with NaN
    matrix = matrix.dropna()
    #
    # df_on = df_clean.pivot_table(index=index,
    #                             columns='uv_components', 
    #                             values='projection_score',
    #                             aggfunc='mean').reset_index()
    #
    # df_off = df_clean.pivot_table(index=index,
    #                             columns='uv_components', 
    #                             values='perpendicular_score',
    #                             aggfunc='mean').reset_index()
    #
    # df_on_abs = df_clean.pivot_table(index=index,
    #                             columns='uv_components', 
    #                             values='abs_projection_score',
    #                             aggfunc='mean').reset_index()
    # df_off_a1bs = df_clean.pivot_table(index=index,
    #                             columns='uv_components', 
    #                             values='abs_perpendicular_score',
    #                             aggfunc='mean').reset_index()
    #
    # import pdb; pdb.set_trace()
    #
    # df_matrix = df_matrix.merge(df_on, on=index)
    # df_matrix = df_matrix.merge(df_off, on=index)
    # df_matrix = df_matrix.merge(df_on_abs, on=index)
    # df_matrix = df_matrix.merge(df_off_a1bs, on=index)
    # df_matrix = df_matrix.dropna()

    df_matrix.set_index(index, inplace=True)
    return df_matrix

def cca_align_animals(df_matrix, index):
    from sklearn.cross_decomposition import CCA
    from sklearn.preprocessing import StandardScaler
    # ------------ CCA ------------
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
    return df_matrix_projected

def prep_uv_magnitude_over_time(df_matrix, time='events'):
    """
    """
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
    def assign_time(group, time=time):
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
    # Define window size
    window_size = 20
    # Apply rolling mean
    df_matrix['magnitude_u_smooth'] = df_matrix.groupby(['animal', 'genH'])['magnitude_u'].transform(lambda x: x.rolling(window_size).mean())
    df_matrix['magnitude_v_smooth'] = df_matrix.groupby(['animal', 'genH'])['magnitude_v'].transform(lambda x: x.rolling(window_size).mean())
    df_matrix_time.head()

    # Subtract minimum time per animal
    df_matrix = df_matrix_time.reset_index()
    df_matrix['time'] = df_matrix['time'] - df_matrix.groupby('animal')['time'].transform('min')

    # Define number of bins
    num_bins = 100
    # Bin the event_time data
    df_matrix['time_bin'] = df_matrix.groupby(['animal', 'genH', 'epoch'])['time'].transform(lambda x: pd.cut(x, num_bins, labels=range(num_bins)))
    return df_matrix


# Load the provided CSV file
df = pd.read_parquet(origin)
if zscore:
    df.query('zscore == True', inplace=True)
else:
    df.query('zscore == False', inplace=True)
# Print the modificiation time of the origin file as XX-XX-XXXX XX:XX:XX
import arrow
print("Origin file last modified on {}".format(arrow.get(os.path.getmtime(origin)).format('MM-DD-YYYY HH:mm:ss')))
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
def commsub_project(df):
    u = df['event_u_values']
    v = df['event_v_values']
    return (u + v) / np.sqrt(2)
def label_epochs(df, time_col='event_time', time_threshold=60*10):
    # Order
    df = df.sort_values(time_col)
    # Calculate time difference
    df.loc[:,'time_diff'] = df[time_col].diff()
    # Initialize the epoch column with 0
    df.loc[:,'epoch2'] = 0
    # If time difference is more than threshold, increment epoch
    df.loc[df['time_diff'] > time_threshold, 'epoch'] = 1
    # Calculate the cumulative sum of the epoch column
    df.loc[:, 'epoch2'] = df['epoch2'].cumsum()
    # Drop the time_diff column as we don't need it anymore
    df.drop(columns=['time_diff'], inplace=True)
    return df
df = (df.sort_values(['event_time','animal'])
      .groupby('animal').progress_apply(label_epochs)
      .rename(columns={'animal': 'animal_new'})
      .reset_index())

# Remove rows with NaN in 'event_u_values' or 'event_v_values'
# df_clean = df.dropna(subset=['event_u_values', 'event_v_values'])
df_clean = df.copy()
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
# df_clean['on_commsub'] = 1 - df_clean['distance_to_line']
df_clean['on_commsub'] = commsub_project(df)
df_clean['on_commsub_mag'] = df_clean['on_commsub'] * df_clean['distance_to_origin']
df_summary = df_clean.groupby(['genH', 'patterns', 'uv_components', 'animal', 'events']).mean().reset_index()
df_clean['off_commsub_mag'] = df_clean['distance_to_line'].abs()
df_clean['on_over_off'] = df_clean['on_commsub_mag'] / df_clean['distance_to_origin'].abs()
if df.patterns.max() == 9:
    mapping = {1:'high', 2:'high', 3:'high', 4:'low', 5:'low', 6:'low', 7:'mid', 8:'mid', 9:'mid'}
    df_clean['highlow'] = df_clean['patterns'].map(mapping)
else:
    mapping = {1:'high', 2:'high', 3:'high', 4:'low', 5:'low', 6:'low'}
    df_clean['highlow'] = df_clean['patterns'].map(mapping)
# Create a composite column combining 'genH' and 'highlow'
df_clean['genH_highlow'] = df_clean['genH'].astype(str) + "_" + df_clean['highlow']
df_clean['highlow_genH'] = df_clean['highlow'] + "_" + df_clean['genH'].astype(str)
df_clean['magnitude_u'] = np.abs(df_clean['event_u_values'])
df_clean['magnitude_v'] = np.abs(df_clean['event_v_values'])
df_clean[['magnitude_u', 'magnitude_v']].head()

overall_pattern = df_clean.patterns.max() + 1
mapping = {int(i):(int(i)-1)%3+1 for i in range(1, int(df.patterns.max()+1))}
mapping[overall_pattern] = np.nan
df_clean['pattern_class'] = df_clean['patterns'].map(mapping)
mapping = {1:'theta', 2:'delta', 3:'ripple'}
df_clean['pattern_class_name'] = df_clean['pattern_class'].map(mapping)

# FILTER # WARNING: if you use 'match' below, you must comment this out
# df_clean = df_clean.query(f'pattern_cca1 == 2 & pattern_cca2 == {df.patterns.max()}')

# Display the first few rows of the dataframe
df_clean.head()
uv_components = sorted(df_clean['uv_components'].unique())
patterns      = sorted(df_clean['patterns'].unique())
genH_values   = df_clean['genH'].unique()
genH_highlow  = sorted(df_clean['genH_highlow'].unique())
print(df_clean.uv_components.unique())


# Calculate the projection score onto the x=y line
df_clean['projection_score'] = (df_clean['event_u_values'] + df_clean['event_v_values']) / np.sqrt(2)
# Calculate the score for distance perpendicular to the x=y line
df_clean['perpendicular_score'] = (df_clean['event_u_values'] - df_clean['event_v_values']) / np.sqrt(2)
df_clean[['projection_score', 'perpendicular_score']].head()
df_clean['abs_projection_score'] = df_clean['projection_score'].abs()
df_clean['abs_perpendicular_score'] = df_clean['perpendicular_score'].abs()
df_clean['proj_over_perp'] = df_clean['abs_projection_score'] / df_clean['abs_perpendicular_score']
df_clean['abs_proj_over_perp'] = df_clean['abs_projection_score'] / df_clean['abs_perpendicular_score']
index = ['animal','genH', 'genH_highlow', 'highlow', 'patterns', 'event_time']
descfield = lambda x: df_clean.groupby(['genH', 'highlow'])[x].describe()


# BUG: Subset
df_clean.query('pattern_cca1 == 2', inplace=True) # 1 is hpchpc
if match_mode:
    tmp = []
    for i in range(1, int(df_clean.patterns.max())):
        tmp.append(df_clean.query(f'pattern_cca2 == {i} & patterns == {i}').copy())
    df_match = pd.concat(tmp)
    df_match.loc[:,"pattern_cca2"] = 0
    del tmp
    df_clean = pd.concat(( df_clean, df_match ), axis=0, ignore_index=True)
    del df_match
    df_clean.query(f"pattern_cca2 == 0 | pattern_cca2 == {df_clean.pattern_cca2.max()}", inplace=True)
else:
    df_clean = df_clean.query(f'pattern_cca2 == {df_clean.pattern_cca2.max()}')


df_matrix = prep_uv_melt(df_clean)
# df_matrix_projected = prep_uv_melt(df_clean, project=True)
df_matrix.dropna(inplace=True)
for i in range(1, 6):
    i = float(i)
    df_matrix[f"r_{i}"] = (df_matrix[f"{i}_u"] + df_matrix[f"{i}_v"]) / np.sqrt(2)
    df_matrix[f"p_{i}"] = (df_matrix[f"{i}_u"] - df_matrix[f"{i}_v"]) / np.sqrt(2)
    df_matrix[f"r_{i}_abs"] = df_matrix[f"r_{i}"].abs()
    df_matrix[f"p_{i}_abs"] = df_matrix[f"p_{i}"].abs()
    df_matrix[f"r_over_p_{i}"] = df_matrix[f"r_{i}_abs"] / df_matrix[f"p_{i}_abs"]
    df_matrix[f"u_{i}"] = df_matrix[f"{i}_u"]
    df_matrix[f"v_{i}"] = df_matrix[f"{i}_v"]


# Let's create histograms of the projection score and perpendicular score
# for genH_name and for highlow
tmp = df_clean.query("highlow != 'mid' & pattern_class == 1").sample(20_000)
g=sns.FacetGrid(data=tmp, col='genH', hue='highlow', height=4, aspect=1)
g.map(sns.histplot, 'abs_projection_score', alpha=0.5, stat='density', cumulative=True)
g.map(lambda *pos, **kws: plt.axvline(0, color='black', linestyle='--'))
g.add_legend()
for ax in g.axes.flatten():
    ax.set_ylim(0, 1)
    ax.set_xlim(-10, 10)
g=sns.FacetGrid(data=tmp, col='genH', hue='highlow', height=4, aspect=1)
g.map(sns.histplot, 'abs_perpendicular_score', alpha=0.5, stat='density', cumulative=True)
g.map(lambda *pos, **kws: plt.axvline(0, color='black', linestyle='--'))
g.add_legend()
for ax in g.axes.flatten():
    ax.set_ylim(0, 1)
    ax.set_xlim(-5, 5)

df_clean.groupby(['genH', 'highlow'])['projection_score'].apply(lambda x: x.skew())
df_clean.groupby(['genH', 'highlow'])['perpendicular_score'].apply(lambda x: x.skew())
df_clean.groupby(['genH', 'highlow'])['projection_score'].apply(lambda x: x.kurtosis())
df_clean.groupby(['genH', 'highlow'])['perpendicular_score'].apply(lambda x: x.kurtosis())



# Import the ipython run magic as a function
import glob
from IPython import get_ipython
ipython = get_ipython()
# Define a function to run a cell in the notebook
ipython.magic("load_ext autoreload")
ipython.magic("autoreload 2")
pyfiles = glob.glob("/Volumes/MATLAB-Drive/Shared/Notebooks/python/eventuv_plots/*.py")
for pyfile in pyfiles:
    ipython.magic("run -i " + pyfile)
