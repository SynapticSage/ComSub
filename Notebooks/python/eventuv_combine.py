import pandas as pd
import os
from glob import glob
import arrow
import datetime

# Directory where your files are located
intermediate = "midpattern=true"
directory = f"/Volumes/MATLAB-Drive/Shared/figures/{intermediate}/tables/"

# List of file names
files = glob(os.path.join(directory, "eventuv*_*.parquet"))

# Empty list to hold dataframes
dfs = []

# Read each file into a dataframe and append to list
for file in files:
    print("Reading {}".format(file))
    ctime = os.path.getctime(file)
    # Convert ctime to human readable format
    ctime = datetime.datetime.fromtimestamp(ctime).strftime('%Y-%m-%d %H:%M:%S')
    print("...created on {}".format(ctime))
    path = os.path.join(directory, file)
    df = pd.read_parquet(path)
    dfs.append(df)

# Combine all dataframes into one
combined_df = pd.concat(dfs, ignore_index=True)

def label_epochs(df, time_col='time', time_threshold=4):
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
if "time" in combined_df.columns:
    combined_df = combined_df.sort_values(['time','animal']).groupby('animal').apply(label_epochs).reset_index()

# Save the combined dataframe as a parquet file
combined_df.to_parquet(os.path.join(directory, "eventuv.parquet"))
combined_df.to_pickle(os.path.join(directory, "eventuv.pickle"))

# WARNINGS:
