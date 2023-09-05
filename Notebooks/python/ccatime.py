import numpy as np
from pathos.helpers import cpu_count
from jax import numpy as jnp
import cupy as cp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from statsmodels.tsa.stattools import grangercausalitytests
from scipy.stats import f
from sklearn.model_selection import TimeSeriesSplit, KFold
import time, os
from PyIF import te_compute as te
from tqdm import tqdm
tqdm.pandas()
from scipy.stats import f
import multiprocessing, dill
intermediate = "midpattern=true"
folder = f'/Volumes/MATLAB-Drive/Shared/figures/{intermediate}/tables/'
name   = 'powerccatime'
def directionality_test(granger_func, x, y, lag_order, xtest=None, ytest=None):
    f_stat_xy, p_value_xy, r2_xy, r2_y, r2_xy_p = granger_func(x, y, lag_order, xtest, ytest)
    f_stat_yx, p_value_yx, r2_yx, r2_x, r2_yx_p = granger_func(y, x, lag_order, xtest, ytest)
    direction = "x->y" if f_stat_xy > f_stat_yx else "y->x"
    significance = "significant" if min(p_value_xy, p_value_yx) < 0.05 else "not significant"
    magnitude = abs(r2_xy - r2_yx)
    magnitude_p = abs(r2_xy_p - r2_yx_p)
    return direction, significance, magnitude, magnitude_p
def is_continuos_var(df, col):
    return df[col].nunique() > 20 
def get_cont_vars(df):
    return [col for col in df.columns if is_continuos_var(df, col)]

# # Example usage
# x = np.random.randn(100)
# y = 0.5 * x + np.random.randn(100)
# lag_order = 5
# start = time.time()
# result  = granger_causality_cupy(x, y, lag_order)
# end = time.time()
# result2,p2 = granger_causality_jax(x, y, lag_order)
# end2 = time.time()
# result3 = grangercausalitytests(np.column_stack([x,y]), lag_order, verbose=False)
# end3 = time.time()
# print(result)
# print(result2)
# print(result3)
# print("cupy")
# print("time: ", end-start)
# print("jax")
# print("time: ", end2-end)
# print("statsmodels")
# print("time: ", end3-end2)

# Load the data
print("Loading data...")
df = pd.read_parquet(os.path.join(folder, f'{name}.parquet'))
df = df[df.nameFull.isnull()]
print("Data loaded.")
cont_vars = get_cont_vars(df)
dfc = df[["animal", *cont_vars]]
print("Drop[ping time...")
cont_vars.remove("time")
dfc = dfc.drop(["time"],axis=1)
# add these columns back from df: 'rewarded', 'trajbound', 'inBoundChoiceTimes', 'outBoundChoiceTimes', 'rewardTimes'
# dfc = pd.concat([dfc, df[['rewarded', 'trajbound', 'inBoundChoiceTimes', 'outBoundChoiceTimes', 'rewardTimes','animal']]], axis=1)
# zscore dfc columns
print("Z-scoring...")
dfc.loc[:,cont_vars] = dfc.loc[:,cont_vars].apply(lambda x: (x - x.mean()) / x.std())
print("Z-scoring done.")
# Define the lag
maxlag = 20
animals = df['animal'].unique()
chunk_size = 260_000 // cpu_count()
i_list_startswith = ["U","V","S","Cavg","wpli_avg"]

def process_pair(pair, k=1, max_lag=maxlag, safetyCheck=False, GPU=False, 
                 calc_te=True, chunk_size=chunk_size):
    i, j = pair
    if i == j or not any(dfc.columns[i].startswith(x) for x in i_list_startswith):
        return None
    print(f"Processing {i}, {j}")
    if dfc.columns[i] == "animal" or dfc.columns[j] == "animal":
        return None
    if any(dfc.iloc[:,[i,j]].dtypes == "object"):
        return None
    RESULTS = []
    for animal in tqdm(dfc['animal'].unique(),total=len(animals)):  # loop through each unique animal
        dfc_animal = dfc.query("animal == @animal")  # get the data for the current animal
        X_animal = dfc_animal.iloc[:, [i, j]].dropna().values
        num_chunks = X_animal.shape[0] // chunk_size  # calculate num_chunks for the current animal
        for chunk_index in range(num_chunks): 
            start_index = chunk_index * chunk_size
            end_index = min(start_index + chunk_size, len(X_animal))
            X_chunk = X_animal[start_index:end_index]
            try:
                results = grangercausalitytests(X_chunk, maxlag=max_lag, verbose=False)
            except ValueError:
                print(f"ValueError for {i}, {j}, {animal}, {chunk_index}")
                continue
            df_results = pd.DataFrame({
                'lag': range(1, max_lag+1),
                'pvalue': [results[lag][0]['ssr_ftest'][1] for lag in range(1, max_lag+1)],
                'F': [results[lag][0]['ssr_ftest'][0] for lag in range(1, max_lag+1)],
                'ssr_chi2test': [results[lag][0]['ssr_chi2test'][0] for lag in range(1, max_lag+1)],
                'ssr_chi2test_pvalue': [results[lag][0]['ssr_chi2test'][1] for lag in range(1, max_lag+1)],
                'lrtest': [results[lag][0]['lrtest'][0] for lag in range(1, max_lag+1)],
                'lrtest_pvalue': [results[lag][0]['lrtest'][1] for lag in range(1, max_lag+1)],
                'params_ftest': [results[lag][0]['params_ftest'][0] for lag in range(1, max_lag+1)],
                'params_ftest_pvalue': [results[lag][0]['params_ftest'][1] for lag in range(1, max_lag+1)],
                'animal': animal,  # add 'animal' to the results
            })
            if calc_te:
                df_results['transfer_entropy'] = te.te_compute(X_chunk[:, 0], X_chunk[:, 1], k, max_lag, safetyCheck, GPU)
                # Compute surrogate transfer entropy with negative lag
                df_results['surrogate_transfer_entropy'] = te.te_compute(X_chunk[::-1, 0], X_chunk[::-1, 1], k, max_lag, safetyCheck, GPU)
            df_results['column1'] = dfc.columns[i]
            df_results['column2'] = dfc.columns[j]
            df_results['chunk_index'] = chunk_index
            RESULTS.append(df_results)
    RESULTS = pd.concat(RESULTS)
    return RESULTS

import gc
gc.collect()
test = process_pair((1,2))
# serialized_func = dill.dumps(process_pair)

## BLOCK DESIGN
# Prepare a list of all pairs of column indices
from pathos.multiprocessing import ProcessingPool as Pool
# Define block size
def process_block(block_pairs):
    with Pool(cpu_count()) as pool:
        block_results = pool.map(process_pair, block_pairs)
    return block_results
pairs = [(i, j) for i in range(dfc.shape[1]) for j in range(dfc.shape[1]) if i != j]

block_size = 32
num_blocks = len(pairs) // block_size + (len(pairs) % block_size != 0)
all_results = []

for block_index in tqdm(range(num_blocks), total=num_blocks, desc='Blocks'):
    if block_index <= 0:
        continue
    start_index = block_index * block_size
    end_index = min((block_index + 1) * block_size, len(pairs))
    block_pairs = pairs[start_index:end_index]
    block_results = process_block(block_pairs)
    all_results.extend(block_results)
    print(f'Finished processing block {block_index + 1} of {num_blocks}')

DF_results = pd.concat(all_results)
pre = name.replace('ccatime', f'_granger_causality')
DF_results.to_csv(os.path.join(folder, f'{pre}.csv'), index=False)

######################################

def find_missing_pairs(DF_results, pairs):
    # Extract indices of processed pairs
    processed_pairs = set(tuple(pair) for pair in DF_results[['column1', 'column2']].drop_duplicates().apply(lambda x: (dfc.columns.get_loc(x['column1']), dfc.columns.get_loc(x['column2'])), axis=1).values)
    
    original_pairs = set(pairs)
    
    missing_pairs = original_pairs - processed_pairs
    
    return missing_pairs

missing_pairs = find_missing_pairs(DF_results, pairs)
print(f"Number of missing pairs: {len(missing_pairs)}")
if missing_pairs:
    print("Missing pairs are:")
    for pair in missing_pairs:
        print(pair)

print("Fraction of missing pairs: ", len(missing_pairs) / len(pairs))

block_size = 32
num_blocks = len(missing_pairs) // block_size + (len(missing_pairs) % block_size != 0)
all_results = []

for block_index in tqdm(range(num_blocks), total=num_blocks, desc='Blocks'):
    if block_index <= 0:
        continue
    start_index = block_index * block_size
    end_index = min((block_index + 1) * block_size, len(missing_pairs))
    block_pairs = list(missing_pairs)[start_index:end_index]
    block_results = process_block(block_pairs)
    all_results.extend(block_results)
    print(f'Finished processing block {block_index + 1} of {num_blocks}')

tmp = pd.concat(all_results)
DF_results = pd.concat([DF_results, tmp])
pre = name.replace('ccatime', f'_granger_causality')
DF_results.to_csv(os.path.join(folder, f'{pre}.csv'), index=False)


######################################

# Let's add some summary statistics
pre = name.replace('ccatime', f'_granger_causality')
DF_results = pd.read_csv(os.path.join(folder, f'{pre}.csv'))
df = DF_results
df['uid'] = (df['chunk_index'].astype(str) + '_' + df['lag'].astype(str) + '_'
             + df['column1'] + '_' +
             df['column2'])
row_dict = df.set_index('uid').T.to_dict()
def directionality_test_vectorized(df, row_dict, drop_missing_reverse=False):
    print("Running directionality test")
    
    uid_reverse = (df['chunk_index'].astype(str) + '_' + df['lag'].astype(str)
                   + '_' + df['column2'] + '_'
                   + df['column1'])
    
    # Get the reverse dataframe
    reverse_rows = [row_dict.get(uid, None) for uid in uid_reverse]
    
    # Check for missing reverse entries
    missing_mask = [row is None for row in reverse_rows]
    if drop_missing_reverse:
        df = df[~np.array(missing_mask)].reset_index(drop=True)
        reverse_rows = [row for row in reverse_rows if row is not None]

    df_reverse = pd.DataFrame(reverse_rows)
    
    f_stat_xy, p_value_xy = (df['F'], df['pvalue'])
    f_stat_yx, p_value_yx = (df_reverse['F'], df_reverse['pvalue'])
    
    direction = np.where(f_stat_xy > f_stat_yx, "x->y", "y->x")
    significance = np.where(np.minimum(p_value_xy, p_value_yx) < 0.05,
                            "significant", "not significant")
    magnitude = np.abs(f_stat_xy - f_stat_yx)
    
    return pd.DataFrame({'direction': direction, 'significance': significance, 'magnitude': magnitude})

df_results = directionality_test_vectorized(DF_results, row_dict, drop_missing_reverse=True)
assert DF_results.shape[0] == df.shape[0]
DF_results = pd.concat([DF_results, df_results], axis=1)
print(f"saving {pre}_directionality.csv")
DF_results.to_csv(os.path.join(folder, f'{pre}_directionality.csv'),
                  index=False)

#  . .     . . .o         |          o          
# -+-+-    | | |.,---.,---|,---.. . ..,---.,---.
# -+-+-    | | |||   ||   ||   || | |||   ||   |
#  ` `     `-'-'``   '`---'`---'`-'-'``   '`---|
#                                          `---'

def capture_windows(df, trigger_column, quantile, w):

    # Remove df columns that are not numeric
    df = df.select_dtypes(include=[np.number])

    # Normalize all columns to be between 0 and 1
    df_normalized = (df - df.min()) / (df.max() - df.min())

    # Compute the threshold
    threshold = df[trigger_column].quantile(quantile)

    # Identify the rows where the trigger column crosses the threshold
    crossing_indices = df.index[df[trigger_column] > threshold]

    # Capture the windows
    windows = []
    for i in crossing_indices:
        if i-w >= 0 and i+w < len(df):
            # Check for gaps greater than 0.2 seconds in the time diff
            window = df_normalized.iloc[i-w:i+w]
            time_diff = window['time'].diff().abs()
            if not any(time_diff > 0.2):
                windows.append(window)

    return windows

wins = capture_windows(df, 'U1', 0.9, 15)


def plot_windows_heatmap(windows):
    # Compute the mean of each column in each window
    window_means = pd.DataFrame([window.mean() for window in windows])

    # Plot the mean of each column over time
    sns.heatmap(window_means.T, vmin=0, center=0.5 vmax=1, cmap='coolwarm')

def plot_windows(windows, quantile, w):
    fig, axs = plt.subplots(len(windows.columns), 1, figsize=(10, 20), constrained_layout=True)
    fig.suptitle(f'Windows for quantile {quantile} and window size {w}', fontsize=16)
    
    for i, col in enumerate(windows.columns):
        data = windows[col]
        mean = data.mean(axis=0)
        # compute 95% confidence interval (1.96 is approximately the 97.5 percentile of a standard normal distribution)
        ci = 1.96 * data.std(axis=0) / np.sqrt(data.shape[0])  
        axs[i].plot(mean, label=f'Mean {col}')
        axs[i].fill_between(range(len(mean)), (mean-ci), (mean+ci), color='b', alpha=.1)
        axs[i].legend(loc='upper right')
    
    plt.show()

plot_windows(captured_windows, quantile, w)

