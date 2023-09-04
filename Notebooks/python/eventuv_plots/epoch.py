
# -----
# 
# BOOTSTRAPPING match levels
# -----
from tqdm import tqdm

def bootstrap_statistics(df, groups, field, statistic, n_boot=1000, ci=95):
    """
    Compute bootstrapped statistics for the given field grouped by the given groups.
    """
    
    # Placeholder dictionaries for results
    all_bootstraps = {}
    cis = []
    
    # Iterate over each unique group
    G = df.groupby(groups)
    for group_values, group in tqdm(G, total=len(G)):
        
        # Bootstrap
        boot_samples = []
        for _ in range(n_boot):
            sample = group[field].sample(len(group), replace=True)
            
            if statistic == 'mean':
                boot_samples.append(sample.mean())
            elif statistic == 'median':
                boot_samples.append(sample.median())
            # ... can add more statistics as needed
            
        # Store all bootstrapped values
        all_bootstraps[group_values] = boot_samples
        
        # Compute the confidence interval
        mean_stat = np.mean(boot_samples)
        lower_bound = np.percentile(boot_samples, (100-ci)/2)
        upper_bound = np.percentile(boot_samples, 100 - (100-ci)/2)
        
        # Append results
        cis.append(list(group_values) + [mean_stat, lower_bound, upper_bound])
    
        # Convert results to dataframes
    df_stats = pd.DataFrame(all_bootstraps).T
    df_cis = pd.DataFrame(cis, columns=groups + ['mean', 'lower_bound', 'upper_bound'])
    df_cis.set_index(groups, inplace=True)
    df_cis.index.set_names(groups, inplace=True)
    df_stats.index.set_names(groups, inplace=True)
    df_cis.columns.set_names(['statistic'], inplace=True)
    df_stats.columns.set_names(['boot'], inplace=True)
    
    return df_stats, df_cis

# Compute the mean and confidence intervals for each category
def corrected_ci_barplot(boot_stats, groups):
    grouped = boot_stats.stack().reset_index().groupby(groups)
    means = grouped[0].mean()
    lower = grouped[0].apply(lambda x: np.percentile(x, 2.5))
    upper = grouped[0].apply(lambda x: np.percentile(x, 97.5))

    # Create a new dataframe with the computed statistics
    stats_df = pd.DataFrame({
        'mean': means,
        'lower': lower,
        'upper': upper
    }).reset_index()

    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111)

    # Define the width of the bars
    width = 0.4

    # Unique values of 'epoch' and 'genH'
    epochs = stats_df['epoch'].unique()
    genH_vals = stats_df['genH'].unique()

    for i, epoch in enumerate(epochs):
        for j, genH_val in enumerate(genH_vals):
            # Filter the stats_df for the current epoch and genH_val
            data_row = stats_df[(stats_df['epoch'] == epoch) & (stats_df['genH'] == genH_val)]
            if not data_row.empty:
                # Calculate position for the bar
                position = i - width/2 + j*width
                # Plot the bar
                ax.bar(position, data_row['mean'].values[0], width=width, label=f"{genH_val if i == 0 else str()}", color=sns.color_palette()[j])
                # Add the error bar
                ax.errorbar(position, data_row['mean'].values[0], yerr=[[data_row['mean'].values[0] - data_row['lower'].values[0]], [data_row['upper'].values[0] - data_row['mean'].values[0]]], fmt='none', color='black')

    # Set the x-ticks and labels
    ax.set_xticks(np.arange(len(epochs)))
    ax.set_xticklabels(epochs)
    ax.set_xlabel('Epoch')
    ax.set_ylabel('Mean Value')
    ax.legend(title='genH')
    ax.title.set_text('Mean Value with 95% Confidence Interval')
    plt.tight_layout()
    plt.show()

def corrected_ci_barplot_v3(df, groups, hue, row=None, col=None, field=None):
    grouped = df.stack().reset_index().groupby(groups)
    means = grouped[0].mean()
    lower = grouped[0].apply(lambda x: np.percentile(x, 2.5))
    upper = grouped[0].apply(lambda x: np.percentile(x, 97.5))

    # Create a new dataframe with the computed statistics
    stats_df = pd.DataFrame({
        'mean': means,
        'lower': lower,
        'upper': upper
    }).reset_index()

    # Define the width of the bars
    width = 0.4

    # Unique values of 'epoch', 'hue', 'row', and 'col'
    epochs = stats_df['epoch'].unique()
    hue_vals = stats_df[hue].unique()
    row_vals = [None] if row is None else stats_df[row].unique()
    col_vals = [None] if col is None else stats_df[col].unique()

    # Create subplots based on the unique values of 'row' and 'col'
    fig, axes = plt.subplots(nrows=len(row_vals), ncols=len(col_vals), figsize=(10*len(col_vals), 5*len(row_vals)))

    # Ensure axes is always a 2D array
    if row is None and col is None:
        axes = np.array([[axes]])
    elif row is None:
        axes = np.array([axes])
    elif col is None:
        axes = np.array([axes]).T

    for r, row_val in enumerate(row_vals):
        for c, col_val in enumerate(col_vals):
            ax = axes[r, c]
            
            for i, epoch in enumerate(epochs):
                for j, hue_val in enumerate(hue_vals):
                    # Filter for the current combination of epoch, hue, row, and col
                    mask = (stats_df['epoch'] == epoch) & (stats_df[hue] == hue_val)
                    if row:
                        mask &= (stats_df[row] == row_val)
                    if col:
                        mask &= (stats_df[col] == col_val)
                    data_row = stats_df[mask]
                    
                    # Calculate position for the bar
                    position = i - width/2 + j*width

                    # Plot the bar
                    if data_row.empty:
                        ax.bar(position, 0, width=width, 
                               label=f"{hue_val if i == 0 else ''}", color=sns.color_palette()[j])
                    else:
                        ax.bar(position, data_row['mean'].values[0], width=width, 
                               label=f"{hue_val if i == 0 else ''}", color=sns.color_palette()[j])
                    
                    # Add the error bar
                    if data_row.empty:
                        ax.errorbar(position, 0, 
                                    yerr=[[0], [0]], 
                                    fmt='none', color='black')
                    else:
                        lower=data_row['mean'].values[0] - data_row['lower'].values[0]
                        upper=data_row['upper'].values[0] - data_row['mean'].values[0]
                        if lower < 0:
                            continue
                        if upper < 0:
                            continue
                        ax.errorbar(position, data_row['mean'].values[0], 
                                    yerr=[[lower], 
                                          [upper]], 
                                    fmt='none', color='black')
            
            # Set the x-ticks, labels, and title for the subplot
            ax.set_xticks(np.arange(len(epochs)))
            ax.set_xticklabels(epochs)
            ax.set_xlabel('Epoch')
            ax.set_ylabel('Mean Value')
            title = []
            if row:
                title.append(f'{row} = {row_val}')
            if col:
                title.append(f'{col} = {col_val}')
            ax.set_title(', '.join(title))
            ax.legend(title=hue)

    # Set share y
    for ax in axes.flat:
        ax.get_shared_y_axes().join(*axes.flat)

    plt.tight_layout()
    if field:
        plt.suptitle("Field: " + field)
        plt.subplots_adjust(top=0.9)
    plt.show()

# ry- changed df_matrix to df_clean

df_matrix['magnitude_r'] = np.sqrt(df_matrix['magnitude_u'] ** 2 + df_matrix['magnitude_v'] ** 2)
df_matrix['on_r_commsub'] = (df_matrix['magnitude_u'] + df_matrix['magnitude_v']) / np.sqrt(2)
df_matrix['off_r_commsub'] = (df_matrix['magnitude_u'] - df_matrix['magnitude_v']) / np.sqrt(2)
df_matrix['on_r_commsub_div_total'] = df_matrix['on_r_commsub'] / df_matrix['magnitude_r']
df_matrix['on_versus_off'] = np.log(df_matrix['on_r_commsub'] / df_matrix['off_r_commsub'])
for i in range(1, 6):
    df_matrix[f"magnitude_r_{i}"] = df_matrix[f"{i}.0_u"] ** 2 + df_matrix[f"{i}.0_v"] ** 2
    df_matrix[f"on_r_commsub_{i}"] = (df_matrix[f"{i}.0_u"] + df_matrix[f"{i}.0_v"]) / np.sqrt(2)
    df_matrix[f"off_r_commsub_{i}"] = (df_matrix[f"{i}.0_u"] - df_matrix[f"{i}.0_v"]) / np.sqrt(2)
    df_matrix[f"on_r_commsub_div_total_{i}"] = df_matrix[f"on_r_commsub_{i}"] / df_matrix[f"magnitude_r_{i}"]
    df_matrix[f"on_versus_off_{i}"] = np.log(df_matrix[f"on_r_commsub_{i}"] / df_matrix[f"off_r_commsub_{i}"])


for field in ['magnitude_u', 'magnitude_v', 'magnitude_r', 'on_r_commsub', 'off_r_commsub', 'on_r_commsub_div_total', 'on_versus_off']:
    print(field)
    boot_stats, boot_cis = bootstrap_statistics(df_matrix, ['epoch', 'highlow', 'highlow_genH', 'genH_highlow', 'genH'], field, 'mean', n_boot=1000, ci=95)
    corrected_ci_barplot_v3(boot_stats, ['epoch', 'highlow', 'genH'], 'highlow', 'genH', field=field)
    plt.savefig(os.path.join(figfolder,f'{field}_genH_overtime_hue=highlow_row=genH.png'), dpi=300)
    plt.savefig(os.path.join(figfolder,f'{field}_genH_overtime_hue=highlow_row=genH.pdf'), dpi=300)

for i in range(1, 3):
    for field in [f"{i}.0_u", f"{i}.0_v", f"r_{i}.0", f"p_{i}.0", f"r_over_p_{i}.0", f"r_{i}.0_abs", f"p_{i}.0_abs"]:
        print(field)
        boot_stats, boot_cis = bootstrap_statistics(df_matrix, ['epoch', 'highlow', 'highlow_genH', 'genH_highlow', 'genH'], field, 'mean', n_boot=1000, ci=95)
        corrected_ci_barplot_v3(boot_stats, ['epoch', 'highlow', 'genH'], 'highlow', 'genH', field=field)
        plt.savefig(os.path.join(figfolder,f'{i}_{field}_genH_overtime_hue=highlow_row=genH.png'), dpi=300)
        plt.savefig(os.path.join(figfolder,f'{i}_{field}_genH_overtime_hue=highlow_row=genH.pdf'))

