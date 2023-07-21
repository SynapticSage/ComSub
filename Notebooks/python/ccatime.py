# Analyzes ccatime tabular datasets
import numpy as np
from jax import numpy as jnp
import cupy as cp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from statsmodels.tsa.stattools import grangercausalitytests
from sklearn.model_selection import TimeSeriesSplit, KFold
import time
from tqdm import tqdm
from pyexfiltrator import exfiltrate, exfiltrated
from scipy.stats import combine_pvalues


def granger_causality_jax(x, y, lag_order):
    # Compute the lagged versions of x and y
    X = jnp.column_stack([x[t-lag_order:t] for t in range(lag_order, len(x))])
    Y = y[None, lag_order:]
    X = X.T
    Y = Y.T
    # Estimate the autoregressive models
    model_xy = jnp.linalg.lstsq(X, Y, rcond=None)[0]
    model_y  = jnp.linalg.lstsq(X[:, 1:], Y, rcond=None)[0]
    # Compute the residual sum of squares
    rss_xy = jnp.sum((Y - jnp.dot(X, model_xy))**2)
    rss_y  = jnp.sum((Y - jnp.dot(X[:, 1:], model_y))**2)
    # Compute the F-statistic
    f_stat = ((rss_y - rss_xy) / lag_order) / (rss_xy / (len(Y) - 2 * lag_order))
    # Compute the P-value
    p_value = 1 - jnp.sum(np.random.f(1, lag_order, len(Y) - 2 * lag_order) > f_stat) / len(Y)
    return f_stat), p_value

def granger_causality_cupy(x, y, lag_order, xtest=None, ytest=None):
    x = cp.asarray(x)
    y = cp.asarray(y)
    # Compute the lagged versions of x and y
    X = cp.column_stack([x[t-lag_order:t] for t in range(lag_order, len(x))])
    Y = y[None, lag_order:]
    X = X.T
    Y = Y.T
    # Estimate the autoregressive models
    model_xy = cp.linalg.lstsq(X, Y, rcond=None)[0]
    model_y = cp.linalg.lstsq(X[:, 1:], Y, rcond=None)[0]
    # Compute the residual sum of squares
    if xtest is None:
        yhat_xy = cp.dot(X, model_xy)
        yhat_y  = cp.dot(X[:, 1:], model_y)
    else:
        xtest = cp.asarray(xtest)
        xtest = cp.column_stack([xtest[t-lag_order:t] for t in range(lag_order, len(xtest))])
        xtest = xtest.T
        yhat_xy = cp.dot(xtest, model_xy)
        yhat_y  = cp.dot(xtest[:, 1:], model_y)
    if ytest is None:
        ytest = Y
    else:
        ytest = cp.asarray(ytest)
        ytest = ytest[None, lag_order:]
        ytest = ytest.T
    total   = cp.sum((ytest - cp.mean(ytest))**2)
    rss_xy  = cp.sum((ytest - yhat_xy)**2)
    rss_y   = cp.sum((ytest - yhat_y)**2)
    # Compute the F-statistic
    f_stat = ((rss_y - rss_xy) / lag_order) / (rss_xy / (len(ytest) - 2 * lag_order))
    # Compute the P-value
    p_value = 1 - cp.sum(cp.random.f(1, lag_order, len(ytest) - 2 * lag_order) > f_stat) / len(ytest)
    r2_xy = 1 - rss_xy / total
    r2_y  = 1 - rss_y / total
    # Compute yhat
    return f_stat.get(), p_value.get(), r2_xy.get(), r2_y.get()

def is_continuos_var(df, col):
    return df[col].nunique() > 20 
def get_cont_vars(df):
    return [col for col in df.columns if is_continuos_var(df, col)]

# Example usage
x = np.random.randn(100)
y = 0.5 * x + np.random.randn(100)
lag_order = 2

start = time.time()
result,p1, yhat_xy, yhat_y  = granger_causality_cupy(x, y, lag_order)
end = time.time()
result2,p2 = granger_causality_jax(x, y, lag_order)
end2 = time.time()
print("cupy")
print(result)
print("time: ", end-start)
print("jax")
print(result2)
print("time: ", end2-end)


# Load the data
df = pd.read_csv('/Volumes/MATLAB-Drive/Shared/figures/tables/ZT2coherenceccatime.csv')
cont_vars = get_cont_vars(df)
dfc = df[cont_vars]
# Define the lag
maxlag = 100
method = 'cupy'
gmeth = granger_causality_jax if method == 'jax' else granger_causality_cupy

# Initialize the TimeSeriesSplit object
n_splits = 2
# tscv = TimeSeriesSplit(n_splits=30)
tscv = KFold(n_splits=n_splits, shuffle=False)
# Loop over each pair of columns in the data frame
DF_results = pd.DataFrame()
for i in tqdm(range(dfc.shape[1]), desc='Column', total=df.shape[1]):
    for j in tqdm(range(i+1, dfc.shape[1]), desc='Row'):
        
        # Extract the time series for the pair of columns
        X = dfc.iloc[:, [i, j]].values
        if method == 'cupy':
            cp._default_memory_pool.free_all_blocks()
            X = cp.asarray(X, dtype=cp.float32)
        elif method == 'jax':
            X = jnp.asarray(X)
        
        # Initialize a list to store the results for each split
        results = []
        
        # Perform time series cross-validation
        itest, (train_index, test_index) = next(enumerate(tscv.split(X)))
        # Get last index of train_index
        itest, (train_index, test_index) = list(enumerate(tscv.split(X)))[-1]
        for itest, (train_index, test_index) in tqdm(enumerate(tscv.split(X)),
                                                     desc='Split', 
                                                     total=n_splits):
            
            # Split the data into training and test sets
            X_train, X_test = X[train_index], X[test_index]
            
            if method == 'cupy' or method == 'jax':
                for lag in tqdm(range(1, maxlag+1), desc='Lag', total=maxlag):
                    # Perform the Granger causality test on the training data
                    assert(X_train.shape[1] > 1)
                    assert(X_train.shape[0] > lag)
                    try:
                        tup = gmeth(X_train[:37_000, 0], X_train[:37_000, 1],
                                lag+1, xtest=X_test[:, 0], ytest=X_test[:, 1])
                        # Extract the p-value for the lag of interest and add it to the results list
                        results.append(dict(
                            lag=lag,
                            pvalue=tup[1],
                            F=tup[0],
                            r2_xy=tup[2],
                            r2_y=tup[3],
                            ))
                        print("SUCCESS")
                    except Exception as e:
                        print(e)
                        results.append(None)
                        print("FAIL")
            else:
                # Perform the Granger causality test on the training data
                t1=time.time()
                result = grangercausalitytests(X_train, maxlag=maxlag, verbose=False)
                t2=time.time()
                print("time: ", t2-t1)
            
                # Extract the p-value for the lag of interest and add it to the results list
                pvalue = result[maxlag][0]['ssr_ftest'][1]
                F = result[maxlag][0]['ssr_ftest'][0]
                results.append([list(range(maxlag)), pvalue, F])
        # Calculate the average p-value across all cross-validation splits
        print("EXFILTRATE")
        exfiltrate()
        # Figure out None fraction
        none_frac = np.sum([x is None for x in results]) / len(results)
        print(f"None fraction: {none_frac}")
        # Remove None values
        results = [x for x in results if x is not None]
        df_results = pd.DataFrame(results)
        DF_results.append(df_results)

                                              
#  . .     . . .o         |          o          
# -+-+-    | | |.,---.,---|,---.. . ..,---.,---.
# -+-+-    | | |||   ||   ||   || | |||   ||   |
#  ` `     `-'-'``   '`---'`---'`-'-'``   '`---|
#                                          `---'

def capture_windows(df, trigger_column, quantile, w):
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

