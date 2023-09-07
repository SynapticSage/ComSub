import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import numpy as np
from tqdm import tqdm

# Load the data from the CSV file
file = '/Volumes/MATLAB-Drive/Shared/figures/midpattern=true/tables/fig2_prediction.csv'
# Print when this file was last modified as XX-XX-XXXX XX:XX:XX
print(f'Last modified: {time.ctime(os.path.getmtime(file))}')
# Print creation
print(f'Created: {time.ctime(os.path.getctime(file))}')
data = pd.read_csv(file)
figfolder = os.path.join(*os.path.split(os.path.dirname(file))[0:-1], 'python_predictions')
if not os.path.exists(figfolder):
    os.makedirs(figfolder)
def classify_highlow(name):
    if 'control' in name:
        return 'low'
    elif 'mid' in name:
        return 'mid'
    else:
        return 'high'
def classify_rhythm(name):
    if 'theta' in name:
        return 'theta'
    elif 'delta' in name:
        return 'delta'
    elif 'ripple' in name:
        return 'ripple'
data.loc[:,'highlow'] = data.name.apply(classify_highlow)
data.loc[:,'rhythm'] = data.name.apply(classify_rhythm)
# Display the first few rows of the dataframe
data.head()
pd.set_option('display.max_rows', None)
data.describe().T
data.info()
data.query('highlow != "mid"', inplace=True)
num_samples = 1000
sample_size = 500

# Describing the variables
# - animal : animal name
# - genH : type of network pattern
# - name : frequency of network pattern
# - direction : brain area directionality of network pattern
# - nSource : number of source neurons
# - nTarget : number of target neurons
# - Var7 (iTarget) : index of target neuron
# - Var8 (perf) : prediction performance
# - mu : mean prediction performance
# - dev : std dev of prediction performance
# - iTarg : index of target neuron
# - highlow : high or low performance (or mid)

# ----------
# PREPROCESSING
# ----------

# Filter out rows where 'Var8 (perf)', 'mu', or 'dev' is less than -2
df = data[(data['Var8'] >= -2) & (data['mu'] >= -2) & (data['dev'] >= -2)]
# Display the first few rows of the filtered dataframe
df.loc[:,'genH_highlow'] = df['genH'] + '_' + df['highlow']
df.loc[:,'dir_genH_highlow'] = df['direction'] + '_' + df['genH'] + '_' + df['highlow']
df.loc[:,'dir_genH'] = df['direction'] + '_' + df['genH']
df.head()

# Filter data for 'hpc-pfc' direction
hpc_pfc_data       = df[df['direction'] == 'hpc-pfc']
hp_noncontrol_data = df[(df['direction'] == 'hpc-pfc') & ~df['name'].str.contains('control')]
hp_control_data    = df[(df['direction'] == 'hpc-pfc') & df['name'].str.contains('control')]
hh_noncontrol_data = df[(df['direction'] == 'hpc-hpc') & ~df['name'].str.contains('control')]
hh_control_data    = df[(df['direction'] == 'hpc-hpc') & df['name'].str.contains('control')]
hpc_pfc_data.attrs['name']       = 'hpc-pfc'
hp_noncontrol_data.attrs['name'] = 'hpc-pfc real'
hp_control_data.attrs['name']    = 'hpc-pfc control'
hh_noncontrol_data.attrs['name'] = 'hpc-hpc real'
hh_control_data.attrs['name']    = 'hpc-hpc control'

plt.ion()

# ----------------------------------------
# SHOWS: Unbalanced
# PLOT: Characterize the data
# INFO: HELPFUL
#-----------------------------------------

def characterize(data):
    import matplotlib.cm as cm
    # Set the temporary rcParams values for color sequence
    with plt.rc_context({'axes.prop_cycle': plt.cycler(color=cm.gray(np.linspace(0, 1, 20)))}):
        # Set 
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        # Plot for 'animal' column
        sns.countplot(ax=axes[0, 0], x='animal', data=data)
        axes[0, 0].set_title('Animal')
        # Plot for 'genH' column
        sns.countplot(ax=axes[0, 1], x='genH', data=data)
        axes[0, 1].set_title('Type of Network Pattern (genH)')
        # Plot for 'name' column
        sns.countplot(ax=axes[1, 0], x='name', data=data)
        axes[1, 0].set_title('Frequency of Network Pattern (name)')
        # Plot for 'direction' column
        sns.countplot(ax=axes[1, 1], x='direction', data=data)
        axes[1, 1].set_title('Brain Area Directionality of Network Pattern (direction)')
        plt.tight_layout()
        plt.show()
characterize(df)

plt.savefig(os.path.join(figfolder, 'characterize_data.png'))
plt.savefig(os.path.join(figfolder, 'characterize_data.pdf'))
plt.close()

# WARNING: Animal field requires balancing

# ----------------------------------------

# # Histogram for continuous variables
# fig, axes = plt.subplots(3, 2, figsize=(15, 10))
#
# # Plot for 'nSource' column
# sns.histplot(ax=axes[0, 0], data=data, x='nSource', kde=True)
# axes[0, 0].set_title('Number of Source Neurons (nSource)')
#
# # Plot for 'nTarget' column
# sns.histplot(ax=axes[0, 1], data=data, x='nTarget', kde=True)
# axes[0, 1].set_title('Number of Target Neurons (nTarget)')
#
# # Plot for 'Var7 (iTarget)' column
# sns.histplot(ax=axes[1, 0], data=data, x='Var7', kde=True)
# axes[1, 0].set_title('Index of Target Neuron (iTarget)')
#
# # Plot for 'Var8 (perf)' column
# sns.histplot(ax=axes[1, 1], data=data, x='Var8', kde=True)
# axes[1, 1].set_title('Prediction Performance (perf)')
#
# # Plot for 'mu' column
# sns.histplot(ax=axes[2, 0], data=data, x='mu', kde=True)
# axes[2, 0].set_title('Mean Prediction Performance (mu)')
#
# # Plot for 'dev' column
# sns.histplot(ax=axes[2, 1], data=data, x='dev', kde=True)
# axes[2, 1].set_title('Std Dev of Prediction Performance (dev)')
#
# plt.tight_layout()
# plt.show()
#
# ----------

# Binning continuous variables
def continuous_binned_characterize(data):
    bins = 50

    fig, axes = plt.subplots(3, 2, figsize=(15, 10))

    # Plot for 'nSource' column
    data['nSource'].plot(kind='hist', bins=bins, ax=axes[0, 0])
    axes[0, 0].set_title('Number of Source Neurons (nSource)')

    # Plot for 'nTarget' column
    data['nTarget'].plot(kind='hist', bins=bins, ax=axes[0, 1])
    axes[0, 1].set_title('Number of Target Neurons (nTarget)')

    # Plot for 'Var7 (iTarget)' column
    data['Var7'].plot(kind='hist', bins=bins, ax=axes[1, 0])
    axes[1, 0].set_title('Index of Target Neuron (iTarget)')

    # Plot for 'Var8 (perf)' column
    data['Var8'].plot(kind='hist', bins=bins, ax=axes[1, 1])
    axes[1, 1].set_title('Prediction Performance (perf)')

    # Plot for 'mu' column
    data['mu'].plot(kind='hist', bins=bins, ax=axes[2, 0])
    axes[2, 0].set_title('Mean Prediction Performance (mu)')

    # Plot for 'dev' column
    data['dev'].plot(kind='hist', bins=bins, ax=axes[2, 1])
    axes[2, 1].set_title('Std Dev of Prediction Performance (dev)')

    plt.tight_layout()
    plt.show()

continuous_binned_characterize(hpc_pfc_data)
plt.savefig(os.path.join(figfolder, 'continuous_binned_characterize_hpc_pfc_data.png'))
plt.savefig(os.path.join(figfolder, 'continuous_binned_characterize_hpc_pfc_data.pdf'))
plt.close()


# ----------
# PLOT: Performance Distribution
# INFO: NOT SO HELPFUL

def performance_distribution_split(data, quantile_xlim=0.99):
    # Define the bins
    bins = 100
    # Define the conditions
    conditions = ['direction', ['direction', 'name'], ['direction', 'name', 'genH'], 'genH']
    # Create the plots
    fig, axes = plt.subplots(len(conditions), 1, figsize=(10, 5 * len(conditions)))
    for ax, condition in zip(axes, conditions):
        if isinstance(condition, list):
            title = ', '.join(condition)
            for value in df.groupby(condition).groups:
                subset = df.loc[df.groupby(condition).groups[value]]
                ax.hist(subset['Var8'], bins=bins, alpha=0.5, label=str(value))
        else:
            title = condition
            for value in df[condition].unique():
                subset = df[df[condition] == value]
                ax.hist(subset['Var8'], bins=bins, alpha=0.5, label=str(value))
        
        ax.set_title(f'Performance Distribution Split by {title}')
        ax.set_xlabel('Performance (Var8)')
        ax.set_ylabel('Frequency')
        ax.set_xlim(df['Var8'].quantile(1 - quantile_xlim),
                    df['Var8'].quantile(quantile_xlim))
        # ax.legend()
    plt.tight_layout()
    plt.show()

performance_distribution_split(hpc_pfc_data, quantile_xlim=0.99)
plt.savefig(os.path.join(figfolder, 'performance_distribution_split_data.png'))
plt.savefig(os.path.join(figfolder, 'performance_distribution_split_data.pdf'))

# ----------
# PLOT: Perf dist by genH and direction
# INFO: HELPFUL

def full_dist_boxplot(data):
    # Create a separate plot for performance split by 'genH' and 'direction'
    fig, ax = plt.subplots(figsize=(10, 5))

    # coolwarm palette
    pal = sns.color_palette("coolwarm", 2)
    key = {'low': 0, 'mid': 5, 'high': 11}
    data = data.sort_values(by=['highlow', 'genH', 'direction'], key=lambda x: x.map(key))
    sns.boxplot(ax=ax, x=data['dir_genH'], y=data['Var8'], hue=data['highlow'], palette=pal)

    def calculate_quantiles(data, groups):
        min_25th_quantile = data.groupby(groups)['Var8'].quantile(0.25).min()
        max_75th_quantile = data.groupby(groups)['Var8'].quantile(0.75).max()
        return min_25th_quantile, max_75th_quantile

    min_25th_quantile, max_75th_quantile = calculate_quantiles(data, ['genH', 'direction'])

    ax.set_title('Performance Distribution Split by genH and direction')
    ax.set_xlabel('genH')
    ax.set_ylabel('Performance (Var8)')
    ax.set_ylim(min_25th_quantile, max_75th_quantile)
    ax.legend(title='Direction', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_ylim([-0.1, 0.2])

    plt.tight_layout()
    plt.show()
full_dist_boxplot(df)

plt.savefig(os.path.join(figfolder, 'full_dist_boxplot_hpc_pfc_data.png'))
plt.savefig(os.path.join(figfolder, 'full_dist_boxplot_hpc_pfc_data.pdf'))


# ----------

# Perform Levene's test for homogeneity of variances
levene_test = stats.levene(hpc_pfc_data['Var8'][hpc_pfc_data['genH'] == 'coherence'],
                           hpc_pfc_data['Var8'][hpc_pfc_data['genH'] == 'power'])

# Run ANOVA if Levene's test is passed (p > 0.05)
if levene_test.pvalue > 0.05:
    # Fit the model
    model = ols("Var8 ~ C(genH)", data=hpc_pfc_data).fit()
    # Perform ANOVA
    anova_table = sm.stats.anova_lm(model, typ=2)
else:
    anova_table = None

print(f'Levene\'s test p-value: {levene_test.pvalue}')
print(anova_table)

# ----------

def calculate_means_of_performance(df):
    # Define the number of samples and sample size

    df = hh_noncontrol_data
    name = df.attrs['name']

    # Initialize lists to hold the sample means
    power_means = []
    wpli_means = []
    coherence_means = []

    # Resample and compute means
    for _ in range(num_samples):
        power_sample = df['Var8'][df['genH'] == 'power'].sample(sample_size, replace=True)
        power_means.append(power_sample.mean())
        
        wpli_sample = df['Var8'][df['genH'] == 'wpli'].sample(sample_size, replace=True)
        wpli_means.append(wpli_sample.mean())
        
        coherence_sample = df['Var8'][df['genH'] == 'coherence'].sample(sample_size, replace=True)
        coherence_means.append(coherence_sample.mean())

    # Create a histogram with transparency
    plt.figure(figsize=(10, 5))
    plt.hist(power_means, bins=50, alpha=0.5, label='power')
    plt.hist(wpli_means, bins=50, alpha=0.5, label='wpli')
    plt.hist(coherence_means, bins=50, alpha=0.5, label='coherence')
    plt.title(f'{name}: Distribution of Sample Means for hpc-pfc Direction Split by genH')
    plt.xlabel('Sample Mean of Performance (Var8)')
    plt.ylabel('Frequency')
    plt.legend()
    plt.show()

calculate_means_of_performance(hpc_pfc_data)

plt.savefig(os.path.join(figfolder, 'calculate_means_of_performance_hpc_pfc_data.png'))
plt.savefig(os.path.join(figfolder, 'calculate_means_of_performance_hpc_pfc_data.pdf'))

# ----------

unique_fields = df['dir_genH_highlow'].unique()

# Function to perform bootstrapping
def bootstrap_means(df, column, field, sample_size, num_samples, stat=np.mean):
    means = []
    smallest_num_animal_sample = df.groupby('animal').size().min()
    for _ in range(num_samples):
        query_string = f"dir_genH_highlow == '{field}'"
        sample = (df.groupby('animal').sample(smallest_num_animal_sample))
        sample = (sample
                  .query(query_string)[column]
                  .sample(sample_size, replace=True))
        if stat == np.mean:
            means.append(sample.mean())
        elif stat == np.median:
            means.append(sample.median())
        elif stat == np.std:
            means.append(sample.std())
        elif stat == stats.skew:
            means.append(sample.skew())
        else:
            means.append(stat(sample))
    return means

# Step 2: Initialize subplots
n = len(unique_fields)
fig, axes = plt.subplots(n, n, figsize=(15, 15), sharex=True, sharey=True)

from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def calculate_bootstraps(df, stat=np.mean):
    bootstrap_means_dict = {}
    for field in tqdm(unique_fields):
        bootstrap_means_dict[field] = bootstrap_means(df, 'Var8', field, sample_size, num_samples, stat=stat)
    return bootstrap_means_dict

def plot_bootstrapped_means(bootstrap_means):
    # Pre-compute bootstrap means for each unique field
    if isinstance(bootstrap_means, pd.DataFrame):
        bootstrap_means = calculate_bootstraps(df) 
    elif isinstance(bootstrap_means, dict):
        pass
    unique_fields = list(bootstrap_means.keys())
    # Step 3 and 4: Plot Data and Color Mapping
    cmap = plt.get_cmap('coolwarm', len(unique_fields))
    n = len(unique_fields)

    sort_keys = {'hpc-hpc_power_low': 0,
                 'hpc-hpc_coherence_low': 1,
                 'hpc-pfc_power_low': 2,
                 'hpc-pfc_coherence_low': 3,
                 'hpc-hpc_power_mid': 4,
                 'hpc-hpc_coherence_mid': 5,
                 'hpc-pfc_power_mid': 6,
                 'hpc-pfc_coherence_mid': 7,
                 'hpc-hpc_power_high': 8,
                 'hpc-hpc_coherence_high': 9,
                 'hpc-pfc_power_high': 10,
                 'hpc-pfc_coherence_high': 11
                 }

    unique_fields = sorted(unique_fields, key=lambda x: sort_keys[x])
    fig, axes = plt.subplots(n, n, figsize=(15, 15), sharex=True, sharey=True)
    for i, field_i in tqdm(enumerate(unique_fields), total=n):
        for j, field_j in tqdm(enumerate(unique_fields), total=n):
            ax = axes[i, j]
            if i < j:
                ax.set_visible(False)
                continue

            if i == j:
                # Diagonal subplot: Plot only one histogram in gray
                means_i = bootstrap_means[field_i]
                ax.hist(means_i, bins=50, 
                        #color='gray', 
                        color=mcolors.to_rgba(cmap(i)),
                        alpha=0.5)
            else:
                # Off-diagonal subplot: Plot pair of histograms
                means_i = bootstrap_means[field_i]
                means_j = bootstrap_means[field_j]
                
                ax.hist(means_i, bins=50, color=mcolors.to_rgba(cmap(i)), alpha=0.5, label=field_i)
                ax.hist(means_j, bins=50, 
                        # color='gray', 
                        color=mcolors.to_rgba(cmap(j)),
                        alpha=0.5, label=field_j)

            newline = '\n'
            downarrow = '\u2193'
            uparrow = '\u2191'
            sinewave = '\u223F'
            rep = lambda x : x.replace('_', newline).replace('hpc', 'HPC').replace('pfc', 'PFC').replace('low', downarrow).replace('high', uparrow).replace('power', sinewave).replace('coherence', " ".join((sinewave, sinewave)))
            if i == j:
                ax.set_title(f"{rep(field_j)}")
            if i == n - 1:
                ax.set_xlabel(f"{rep(field_j)}")
            if j == 0:
                ax.set_ylabel(f"{rep(field_i)}")
    # Step 5: Link X-Axes (already done through `sharex=True` and `sharey=True`)
    plt.show()


bootstraps_per_rhythm = {rhythm:calculate_bootstraps(df.query(f'rhythm == "{rhythm}"')) for rhythm in df['rhythm'].unique()}
# let's make this a for loop instead
for rhythm in df['rhythm'].unique():
    plot_bootstrapped_means(bootstraps_per_rhythm[rhythm])
    plt.savefig(os.path.join(figfolder, f'calculate_means_of_performance_{rhythm}_data.png'))
    plt.savefig(os.path.join(figfolder, f'calculate_means_of_performance_{rhythm}_data.pdf'))


# Let's do an upper clipped mean, the mean of the top 10% of the data
def upper_clipped_mean(data, quantile=0.9):
    return data[data > data.quantile(quantile)].mean()
clipped_bootstraps_per_rhythm = {rhythm:calculate_bootstraps(df.query(f'rhythm == "{rhythm}"'), stat=upper_clipped_mean) for rhythm in df['rhythm'].unique()}
# let's make this a for loop instead
for rhythm in df['rhythm'].unique():
    plot_bootstrapped_means(clipped_bootstraps_per_rhythm[rhythm])
    plt.suptitle("Upper Clipped Mean of " + rhythm)
    plt.savefig(os.path.join(figfolder, f'calculate_means_of_performance_{rhythm}_data_clipped.png'))
    plt.savefig(os.path.join(figfolder, f'calculate_means_of_performance_{rhythm}_data_clipped.pdf'))

# Let's do a median
median_bootstraps_per_rhythm = {rhythm:calculate_bootstraps(df.query(f'rhythm == "{rhythm}"'), stat=np.median) for rhythm in df['rhythm'].unique()}
# let's make this a for loop instead
for rhythm in df['rhythm'].unique():
    plot_bootstrapped_means(median_bootstraps_per_rhythm[rhythm])
    plt.suptitle("Median of " + rhythm)
    plt.savefig(os.path.join(figfolder, f'calculate_means_of_performance_{rhythm}_data_median.png'))
    plt.savefig(os.path.join(figfolder, f'calculate_means_of_performance_{rhythm}_data_median.pdf'))


