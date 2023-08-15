import pandas as pd
import os
import glob

# specify the directory you're working in
intermediate = "midpattern=true"
folder = f'/Volumes/MATLAB-Drive/Shared/figures/{intermediate}/tables/'
# find all CSV files in the directory
all_files = glob.glob(os.path.join(folder,'*_faxis=Inf*_regress.csv'))
# filter the files to only those with 'faxis=Inf' in their title
inf_files = [file for file in all_files if 'faxis=Inf' in file]

# create a list to hold the dataframes
df_list = []

# loop through the list of files
for filename in inf_files:
    # read the CSV file into a dataframe
    df = pd.read_csv(filename, index_col=None, header=0)
    
    # create a new column with the filename
    df['filename'] = filename
    
    # append the dataframe to the list
    df_list.append(df)

# concatenate all the dataframes in the list into one dataframe
combined_df = pd.concat(df_list, axis=0, ignore_index=True)

# save the combined dataframe to a new CSV file
combined_df.to_csv(os.path.join(folder,
                                'combined_faxis=Inf_regress.csv'), index=False)
