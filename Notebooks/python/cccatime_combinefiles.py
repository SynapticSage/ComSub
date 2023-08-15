import pandas as pd
import glob
import os
from tqdm import tqdm

# Glob pattern to match the files
intermediate = "midpattern=true"
folder       = f'/Volumes/MATLAB-Drive/Shared/figures/{intermediate}/tables/'
networkpat   = "coherence"
pattern = f"*{networkpat}ccatime.csv"

# Get a list of all matching files
files = glob.glob(os.path.join(folder, pattern))
if len(files) == 0:
    raise ValueError("No files found matching pattern {}".format(pattern))

# Initialize an empty dictionary to hold dataframes
dfs = {}

# Loop over the files and read them into pandas dataframes
for file in tqdm(files, desc="Reading files", total=len(files)):
    # The key will be the rat's name, extracted from the filename
    key = os.path.splitext(file)[0].replace("powerccatime", "")
    key = os.path.basename(key)
    # Read the file into a dataframe and store it in the dictionary
    dfs[key] = pd.read_csv(file)

# Concatenate all the dataframes together
dfs = pd.concat(dfs.values())
dfs.head()
dfs.drop(columns=["animalnames"], inplace=True) # already an animal column

# Save the result to a feather file
dfs.to_parquet(os.path.join(folder, "powerccatime.parquet"))

del dfs
