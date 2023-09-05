import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols

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
data.loc[:,'highlow'] = data.name.apply(lambda x: 'high' if 'control' not in x else 'low')
# Display the first few rows of the dataframe
data.head()

# ----------

# Filter out rows where 'Var8 (perf)', 'mu', or 'dev' is less than -2
filtered_data = data[(data['Var8'] >= -2) & (data['mu'] >= -2) & (data['dev'] >= -2)]
# Display the first few rows of the filtered dataframe
filtered_data.head()
# Filter data for 'hpc-pfc' direction
hpc_pfc_data       = filtered_data[filtered_data['direction'] == 'hpc-pfc']
hp_noncontrol_data = filtered_data[(filtered_data['direction'] == 'hpc-pfc') & ~filtered_data['name'].str.contains('control')]
hp_control_data    = filtered_data[(filtered_data['direction'] == 'hpc-pfc') & filtered_data['name'].str.contains('control')]
hh_noncontrol_data = filtered_data[(filtered_data['direction'] == 'hpc-hpc') & ~filtered_data['name'].str.contains('control')]
hh_control_data    = filtered_data[(filtered_data['direction'] == 'hpc-hpc') & filtered_data['name'].str.contains('control')]
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
characterize(hpc_pfc_data)
plt.savefig(os.path.join(figfolder, 'characterize_hpc_pfc_data.png'))
plt.savefig(os.path.join(figfolder, 'characterize_hpc_pfc_data.pdf'))
plt.close()

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
            for value in filtered_data.groupby(condition).groups:
                subset = filtered_data.loc[filtered_data.groupby(condition).groups[value]]
                ax.hist(subset['Var8'], bins=bins, alpha=0.5, label=str(value))
        else:
            title = condition
            for value in filtered_data[condition].unique():
                subset = filtered_data[filtered_data[condition] == value]
                ax.hist(subset['Var8'], bins=bins, alpha=0.5, label=str(value))
        
        ax.set_title(f'Performance Distribution Split by {title}')
        ax.set_xlabel('Performance (Var8)')
        ax.set_ylabel('Frequency')
        ax.set_xlim(filtered_data['Var8'].quantile(1 - quantile_xlim),
                    filtered_data['Var8'].quantile(quantile_xlim))
        # ax.legend()
    plt.tight_layout()
    plt.show()

performance_distribution_split(hpc_pfc_data, quantile_xlim=0.99)
plt.savefig(os.path.join(figfolder, 'performance_distribution_split_hpc_pfc_data.png'))
plt.savefig(os.path.join(figfolder, 'performance_distribution_split_hpc_pfc_data.pdf'))

# ----------
# PLOT: Perf dist by genH and direction
# INFO: HELPFUL

def full_dist_boxplot(data):
    # Create a separate plot for performance split by 'genH' and 'direction'
    fig, ax = plt.subplots(figsize=(10, 5))

    sns.boxplot(ax=ax, x=filtered_data['genH'], y=filtered_data['Var8'], hue=filtered_data['direction'])

    min_25th_quantile, max_75th_quantile = calculate_quantiles(filtered_data, ['genH', 'direction'])

    ax.set_title('Performance Distribution Split by genH and direction')
    ax.set_xlabel('genH')
    ax.set_ylabel('Performance (Var8)')
    ax.set_ylim(min_25th_quantile, max_75th_quantile)
    ax.legend(title='Direction', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    plt.show()

full_dist_boxplot(hpc_pfc_data)
plt.savefig(os.path.join(figfolder, 'full_dist_boxplot_hpc_pfc_data.png'))
plt.savefig(os.path.join(figfolder, 'full_dist_boxplot_hpc_pfc_data.pdf'))


# ----------

# Perform Levene's test for homogeneity of variances
levene_test = stats.levene(hpc_pfc_data['Var8'][hpc_pfc_data['genH'] == 'coherence'],
                           hpc_pfc_data['Var8'][hpc_pfc_data['genH'] == 'correlation'],
                           hpc_pfc_data['Var8'][hpc_pfc_data['genH'] == 'sparse'])

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

    num_samples = 1000
    sample_size = 500
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

# Store the dataframes in a list
dfs = [hpc_pfc_data, hp_noncontrol_data, hp_control_data, hh_noncontrol_data, hh_control_data]

# Create subplots
fig, axes = plt.subplots(len(dfs), 1, figsize=(10, 5 * len(dfs)))

for df, ax in zip(dfs, axes):
    name = df.attrs['name']

    # Initialize lists to hold the sample means
    power_means = []
    wpli_means = []
    coherence_means = []

    # Resample and compute means
    for _ in range(num_samples):
        power_sample = df['Var8'][df['genH'] == 'power'].sample(sample_size, replace=True)
        power_means.append(power_sample.mean())

        # wpli_sample = df['Var8'][df['genH'] == 'wpli'].sample(sample_size, replace=True)
        # wpli_means.append(wpli_sample.mean())

        coherence_sample = df['Var8'][df['genH'] == 'coherence'].sample(sample_size, replace=True)
        coherence_means.append(coherence_sample.mean())

    # Create a histogram with transparency
    ax.hist(power_means, bins=50, alpha=0.5, label='power')
    # ax.hist(wpli_means, bins=50, alpha=0.5, label='wpli')
    ax.hist(coherence_means, bins=50, alpha=0.5, label='coherence')
    ax.set_title(f'{name}: Distribution of Sample Means')
    ax.set_xlabel('Sample Mean of Performance (Var8)')
    ax.set_ylabel('Frequency')
    ax.legend()

plt.tight_layout()
plt.show()
plt.savefig(os.path.join(figfolder, 'calculate_means_of_performance_all_data.png'))
plt.savefig(os.path.join(figfolder, 'calculate_means_of_performance_all_data.pdf'))
