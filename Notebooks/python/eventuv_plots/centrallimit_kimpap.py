import warnings
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from typing import Callable
import numpy as np
from scipy import stats
from tqdm import tqdm

# default_stat = np.mean
# default_stat = np.std
# default_stat = np.median
# default_stat = stats.skew
# default_stat = stats.kurtosis
def absmean(x):
    return np.mean(np.abs(x))
def absskew(x):
    return stats.skew(np.abs(x))
stat_seq = [
        np.mean, np.std, 
        absmean, 
        np.median, stats.skew, stats.kurtosis
        ]

for default_stat in tqdm(stat_seq, total=6, desc="default_stat"):

    if default_stat == np.mean:
        central_limit_kimpap_folder = os.path.join(figfolder, "central_limit_kimpap")
    else:
        central_limit_kimpap_folder = os.path.join(figfolder, f"central_limit_kimpap_{default_stat.__name__}")
    if not os.path.exists(central_limit_kimpap_folder):
        os.makedirs(central_limit_kimpap_folder)

    # Let's do the same for power high and power low now

    # --------------------


    def bootstrap_means(df_clean: pd.DataFrame, base_query_str: str, additional_queries: list, 
                        statistic: Callable = default_stat, n_samples: int = 100, 
                        simplify_query_info: bool = False) -> pd.DataFrame:
        """
        Bootstrap resampling to estimate the distribution of means for projection_score and perpendicular_score.
        Allows for additional query strings for comparison.
        
        Parameters:
        - df_clean: DataFrame containing the cleaned data
        - base_query_str: Base query string to filter the DataFrame
        - additional_queries: List of additional query strings to be appended to the base query
        - statistic: Function to compute the statistic of interest (default is np.mean)
        - n_samples: Number of bootstrap samples (default is 100)
        - simplify_query_info: Whether to simplify the query info by keeping only the part after '==' (default is False)
        
        Returns:
        - DataFrame containing the bootstrap resampled statistics along with query information
        """
        # Initialize lists to store results
        bootstrap_means_proj = []
        bootstrap_means_perp = []
        query_info = []
        
        # Loop over each additional query string
        for add_query in additional_queries:
            # Combine base query string with additional query string
            full_query_str = f"{base_query_str} & {add_query}"
            
            # Filter the DataFrame based on the full query string
            filtered_df = df_clean.query(full_query_str)
            
            # Perform bootstrap resampling
            for i in range(n_samples):
                # Resample the data
                resampled_df = filtered_df.sample(n=len(filtered_df), replace=True)
                
                # Compute the statistic for projection_score and perpendicular_score
                mean_proj = statistic(resampled_df['projection_score'])
                mean_perp = statistic(resampled_df['perpendicular_score'])
                
                # Store the computed means and query info
                bootstrap_means_proj.append(mean_proj)
                bootstrap_means_perp.append(mean_perp)
                
                # Simplify the query info if specified
                simplified_query = add_query.split('==')[-1].strip() if simplify_query_info else add_query
                simplified_query = simplified_query.replace('"', '') if simplify_query_info else simplified_query
                
                query_info.append(simplified_query)
        
        # Create a DataFrame to store the bootstrap means and query info
        bootstrap_means_df = pd.DataFrame({
            'bootstrap_means_proj': bootstrap_means_proj,
            'bootstrap_means_perp': bootstrap_means_perp,
            'query_info': query_info
        })
        
        return bootstrap_means_df


    def plot_bootstrap_means(bootstrap_df: pd.DataFrame, confidence_interval: float = 0.95, 
                             plot_ci: bool = True, cmap: str = 'viridis', reverse_cmap: bool = False,
                             plot_zero_lines: bool = False, ax=None):
        """
        Plot the bootstrap resampled means using KDE and scatter plots, with optional confidence intervals.
        
        Parameters:
        - bootstrap_df: DataFrame containing the bootstrap resampled means and query_info
        - confidence_interval: Confidence level for confidence intervals (default is 0.95)
        - plot_ci: Whether to plot confidence intervals as axvline and axhline (default is True)
        - cmap: Colormap to use for the plot (default is 'viridis')
        - reverse_cmap: Whether to reverse the colormap (default is False)
        - plot_zero_lines: Whether to plot dark solid black lines at x=0 and y=0 (default is False)
        """
        
        # Unique query_info labels and corresponding colors
        unique_queries = bootstrap_df['query_info'].unique()
        cmap_colors = plt.cm.get_cmap(cmap, len(unique_queries))
        
        if reverse_cmap:
            cmap_colors = cmap_colors.reversed()
        
        # Create a color map for unique query_info labels
        color_map = {query: cmap_colors(i) for i, query in enumerate(unique_queries)}

        badinds = bootstrap_df[['bootstrap_means_proj', 'bootstrap_means_perp']].isna().any(axis=1)
        bootstrap_df = bootstrap_df[~badinds]
        
        if (bootstrap_df.bootstrap_means_proj.unique().shape[0] > 1) and \
           (bootstrap_df.bootstrap_means_perp.unique().shape[0] > 1):
            # KDE plot
            sns.kdeplot(data=bootstrap_df, x="bootstrap_means_proj", y="bootstrap_means_perp", 
                        hue="query_info", fill=True, alpha=0.5, palette=color_map, ax=ax)
            
            # Scatter plot overlay
            sns.scatterplot(data=bootstrap_df, x="bootstrap_means_proj", y="bootstrap_means_perp", 
                            hue="query_info", alpha=0.6, edgecolor=None, palette=color_map, ax=ax)
            
            if plot_ci:
                # Calculate and plot confidence intervals
                for query in unique_queries:
                    subset_df = bootstrap_df[bootstrap_df['query_info'] == query]
                    
                    # Calculate confidence intervals for projection_score and perpendicular_score
                    lower_bound_proj = subset_df['bootstrap_means_proj'].quantile((1 - confidence_interval) / 2)
                    upper_bound_proj = subset_df['bootstrap_means_proj'].quantile(1 - (1 - confidence_interval) / 2)
                    
                    lower_bound_perp = subset_df['bootstrap_means_perp'].quantile((1 - confidence_interval) / 2)
                    upper_bound_perp = subset_df['bootstrap_means_perp'].quantile(1 - (1 - confidence_interval) / 2)
                    
                    # Plot confidence intervals with corresponding colors
                    ax.axvline(x=lower_bound_proj, color=color_map[query], linestyle='--')
                    ax.axvline(x=upper_bound_proj, color=color_map[query], linestyle='--')
                    
                    ax.axhline(y=lower_bound_perp, color=color_map[query], linestyle='--')
                    ax.axhline(y=upper_bound_perp, color=color_map[query], linestyle='--')
        
        # Plot dark solid black lines at x=0 and y=0 if specified
        if plot_zero_lines:
            ax.axvline(x=0, color='black', linestyle='-')
            ax.axhline(y=0, color='black', linestyle='-')
        
        ax.set_xlabel('Bootstrap Means - Projection Score')
        ax.set_ylabel('Bootstrap Means - Perpendicular Score')
        ax.set_title('Bootstrap Resampled Means with Confidence Intervals')
        ax.legend(title='Query Info')
        plt.show()

    # Function to perform a two-sample t-test on raw data from df_clean and return a shorter string with the difference and p-value
    def perform_ttest_raw(df, base_query, query1, query2, column):
        df1 = df.query(f"{base_query} & {query1}")
        df2 = df.query(f"{base_query} & {query2}")
        t_stat, p_val = stats.ttest_ind(df1[column], df2[column])
        return f"Δ {column} = {t_stat:.2f}, p = {p_val:.3f}"
    # Function to perform a two-sample t-test and return a shorter string with the difference and p-value
    def perform_ttest_short(df1, df2, column):
        t_stat, p_val = stats.ttest_ind(df1[column], df2[column])
        return f"Δ {column.split('_')[-1]} = {t_stat:.2f}, p = {p_val:.3f}"


    # ----------------------------------
    def plot_bootstrap_grid(df_clean: pd.DataFrame, base_query_str: str, **kwargs):

        # Generate the bootstrap samples for all the conditions: "coherence_high", "coherence_low", "power_high", "power_low"
        additional_queries_all = ['genH_highlow == "coherence_high"', 'genH_highlow == "coherence_low"', 
                                  'genH_highlow == "power_high"', 'genH_highlow == "power_low"']
        bootstrap_means_df = bootstrap_means(df_clean, base_query_str, additional_queries_all, 
                                                 n_samples=100, simplify_query_info=True)
        bootstrap_means_df.loc[:,'query_info'] =  bootstrap_means_df.query_info.str.replace('"', '')

        # List of unique conditions
        unique_conditions = bootstrap_means_df['query_info'].unique()

        # Initialize the subplot grid
        fig, axes = plt.subplots(4, 4, figsize=(20, 20), sharex=True, sharey=True)
        use_raw_data = False  # Flag to switch between resampled means and raw data

        # Loop through each subplot to generate the plot
        i, row_condition = 0, unique_conditions[0]
        for i, row_condition in tqdm(enumerate(unique_conditions), total=4, desc="row_condition"):
            j, col_condition = 0, unique_conditions[0]
            for j, col_condition in enumerate(unique_conditions):
                ax = axes[i, j]
                
                # Blank out the axes above the diagonal
                if j > i:
                    ax.axis('off')
                    continue
                
                # For off-diagonal plots, use plot_bootstrap_means_in_axes and add t-test results
                if i != j:
                    subset_df1 = bootstrap_means_df.query(f"query_info == '{row_condition}'")
                    subset_df2 = bootstrap_means_df.query(f"query_info == '{col_condition}'")
                    if subset_df1.empty or subset_df2.empty:
                        warnings.warn(f"No data for {row_condition} vs {col_condition}")
                    
                    # Perform t-tests for both axes and get shorter strings
                    if use_raw_data:
                        ttest_proj_str = perform_ttest_raw(df_clean, base_query_str, row_condition, col_condition, 'projection_score')
                        ttest_perp_str = perform_ttest_raw(df_clean, base_query_str, row_condition, col_condition, 'perpendicular_score')
                    else:
                        ttest_proj_str = perform_ttest_short(subset_df1, subset_df2, 'bootstrap_means_proj')
                        ttest_perp_str = perform_ttest_short(subset_df1, subset_df2, 'bootstrap_means_perp')
                    
                    plot_bootstrap_means(subset_df1.append(subset_df2), plot_ci=True, cmap='coolwarm', plot_zero_lines=True, ax=ax, **kwargs)
                    
                    # Add plot titles and labels with t-test results
                    ax.set_title(f"{row_condition} vs {col_condition}\n{ttest_proj_str}\n{ttest_perp_str}")
                
                # For diagonal plots, show the overall KDE and scatter plot for the condition using a constant black colormap
                else:
                    subset_df = bootstrap_means_df[bootstrap_means_df['query_info'] == row_condition]
                    plot_bootstrap_means(subset_df, plot_ci=False, cmap=ListedColormap(['black']), plot_zero_lines=True, ax=ax, **kwargs)
                    
                    # Add plot titles and labels
                    ax.set_title(f"{row_condition}")

        plt.tight_layout()
        plt.show()

    def process_comp(comp, condition=None, N=200_000):

        if condition:
            tmp = df_clean.query(condition)
        else:
            tmp = df_clean
            condition = ""

        if N:
            attempts = 10
            while attempts > 0:
                try:
                    tmp = tmp.groupby(['genH', 'animal']).sample(N)
                    attempts = 0
                except:
                    attempts -= 1
                    print(f"Failed to sample {N} rows. Trying again with {N//2} rows.")
                    N = int(N // 1.5)

        base_query_str =  f'pattern_cca1 == 2 & pattern_cca2==10 & uv_components == {comp}'
        plot_bootstrap_grid(tmp,  base_query_str + " & pattern_class_name == 'theta'")
        plt.suptitle(f"Bootstrap Resampled {default_stat.__name__} with Confidence Intervals - THETA")
        plt.subplots_adjust(top=0.95)
        plt.savefig(os.path.join(central_limit_kimpap_folder, f'plot_bootstrap_grid_theta_uv{comp}{condition}.png'), dpi=300)
        plt.savefig(os.path.join(central_limit_kimpap_folder, f'plot_bootstrap_grid_theta_uv{comp}{condition}.pdf'))

        plot_bootstrap_grid(tmp,  base_query_str + " & pattern_class_name == 'delta'", reverse_cmap=True)
        plt.suptitle(f"Bootstrap Resampled {default_stat.__name__} with Confidence Intervals - delta")
        plt.subplots_adjust(top=0.95)
        plt.savefig(os.path.join(central_limit_kimpap_folder, f'plot_bootstrap_grid_delta_uv{comp}{condition}.png'), dpi=300)
        plt.savefig(os.path.join(central_limit_kimpap_folder, f'plot_bootstrap_grid_delta_uv{comp}{condition}.pdf'))

        plot_bootstrap_grid(tmp,  base_query_str + " & pattern_class_name == 'ripple'")
        plt.suptitle(f"Bootstrap Resampled {default_stat.__name__} with Confidence Intervals - ripple\n" + default_stat.__name__)
        plt.subplots_adjust(top=0.95)
        plt.savefig(os.path.join(central_limit_kimpap_folder, f'plot_bootstrap_grid_ripple_uv{comp}{condition}.png'), dpi=300)
        plt.savefig(os.path.join(central_limit_kimpap_folder, f'plot_bootstrap_grid_ripple_uv{comp}{condition}.pdf'))

    plt.close('all')
    process_comp(1)
    plt.close('all')
    process_comp(2)
    plt.close('all')

    def process_all_comp(condition=None, N=300_000):

        if condition:
            tmp = df_clean.query(condition)
        else:
            tmp = df_clean
            condition = ""

        if N:
            attempts = 10
            while attempts > 0:
                try:
                    tmp = tmp.groupby(['genH', 'animal']).sample(N)
                    attempts = 0
                except:
                    attempts -= 1
                    print(f"Failed to sample {N} rows. Trying again with {N//2} rows.")
                    N = int(N // 2)

        base_query_str =  f'pattern_cca1 == 2 & pattern_cca2==10'
        plot_bootstrap_grid(tmp,  base_query_str + " & pattern_class_name == 'theta'")
        plt.suptitle(f"Bootstrap Resampled {default_stat.__name__} with Confidence Intervals - THETA")
        plt.subplots_adjust(top=0.95)
        plt.savefig(os.path.join(central_limit_kimpap_folder, f'plot_bootstrap_grid_theta_uv{condition}.png'), dpi=300)
        plt.savefig(os.path.join(central_limit_kimpap_folder, f'plot_bootstrap_grid_theta_uv{condition}.pdf'))

        plot_bootstrap_grid(tmp,  base_query_str + " & pattern_class_name == 'delta'", reverse_cmap=True)
        plt.suptitle(f"Bootstrap Resampled {default_stat.__name__} with Confidence Intervals - delta")
        plt.subplots_adjust(top=0.95)
        plt.savefig(os.path.join(central_limit_kimpap_folder, f'plot_bootstrap_grid_delta_uv{condition}.png'), dpi=300)
        plt.savefig(os.path.join(central_limit_kimpap_folder, f'plot_bootstrap_grid_delta_uv{condition}.pdf'))

        plot_bootstrap_grid(tmp,  base_query_str + " & pattern_class_name == 'ripple'")
        plt.suptitle(f"Bootstrap Resampled {default_stat.__name__} with Confidence Intervals - ripple")
        plt.subplots_adjust(top=0.95)
        plt.savefig(os.path.join(central_limit_kimpap_folder, f'plot_bootstrap_grid_ripple_uv{condition}.png'), dpi=300)
        plt.savefig(os.path.join(central_limit_kimpap_folder, f'plot_bootstrap_grid_ripple_uv{condition}.pdf'))

    process_all_comp()

    # - -- 

    def plot_bootstrap_grid_by_epoch(df_clean: pd.DataFrame, base_query_str: str, combine_epochs: bool = True, **kwargs):
        # List of unique conditions
        additional_queries= [
            'genH_highlow == "coherence_high"', 
            'genH_highlow == "coherence_low"', 
            'genH_highlow == "power_high"', 
            'genH_highlow == "power_low"'
        ]

        # Generate the colormap
        unique_epochs  = sort(df_clean['epoch'].dropna().unique())
        cmap = plt.cm.get_cmap('coolwarm')
        n_colors = len(unique_epochs)
        colors = cmap(np.linspace(0, 1, n_colors))

        # Initialize the figure
        fig, axes = plt.subplots(len(additional_queries), 1, figsize=(15, 5 * len(additional_queries)), sharex=True, sharey=True)

        # Loop through each condition
        idx, condition = 0, additional_queries[0]
        for idx, condition in enumerate(additional_queries):
            ax = axes[idx]
            
            # If combine_epochs is True, we combine all epochs into one subplot
            if combine_epochs:
                # Initialize a DataFrame to hold bootstrap_means for all epochs
                combined_bootstrap_means_df = pd.DataFrame()
                
                # Loop through each epoch and get the bootstrap_means
                for epoch_idx, epoch in enumerate(unique_epochs):
                    epoch_query_str = f"{base_query_str} & epoch == {epoch}"
                    bootstrap_means_epoch_df = bootstrap_means(df_clean, epoch_query_str, [condition], 
                                                              n_samples=100, simplify_query_info=True)
                    bootstrap_means_epoch_df['epoch'] = epoch  # Add epoch information to the DataFrame
                    
                    combined_bootstrap_means_df = pd.concat([combined_bootstrap_means_df, bootstrap_means_epoch_df])
                
                # Plot the combined data
                sns.scatterplot(data=combined_bootstrap_means_df, x='bootstrap_means_proj', y='bootstrap_means_perp', 
                                hue='epoch', palette=cmap, ax=ax)
                ax.axvline([0], color='black', linestyle='--')
                ax.axhline([0], color='black', linestyle='--')
            else:
                # Loop through each epoch to create individual plots for each epoch
                for epoch_idx, epoch in enumerate(unique_epochs):
                    epoch_query_str = f"{base_query_str} & epoch == {epoch}"
                    bootstrap_means_epoch_df = bootstrap_means(df_clean, epoch_query_str, [condition], 
                                                              n_samples=100, simplify_query_info=True)
                    sns.scatterplot(data=bootstrap_means_epoch_df, x='bootstrap_means_proj', y='bootstrap_means_perp', 
                                    color=colors[epoch_idx], ax=ax, label=f'Epoch {epoch}')
            ax.axvline([0], color='black', linestyle='--')
            ax.axhline([0], color='black', linestyle='--')
            
            ax.set_title(f"{condition} Over Epochs")
            ax.set_xlabel("Projection Score")
            ax.set_ylabel("Perpendicular Score")
            
        plt.tight_layout()
        plt.show()

    # Renaming the function to better reflect its purpose: Generate and save epoch-based bootstrap plots for different conditions and animals
    def generate_and_save_epoch_bootstrap_plots(df, base_query_str, pattern_class_names, 
                                                animals=None, combine_epochs=True, figappend=""):
        pattern_class = pattern_class_names[0]
        for pattern_class in pattern_class_names:
            animals = animals or ['all']
            for animal in tqdm(animals or ['all'], desc="animal"):
                animal_str = f' & animal == "{animal}"' if animal != 'all' else ''
                query_str = f'{base_query_str} & pattern_class_name == "{pattern_class}"{animal_str}'
                
                plot_bootstrap_grid_by_epoch(df, query_str, combine_epochs)
                
                plt.suptitle(f"Epoch, Bootstrap Resampled {default_stat.__name__} with Confidence Intervals - {pattern_class.upper()} - Animal {animal}")
                
                # Save the plots
                animal_suffix = f"_animal_{animal}" if animal != 'all' else ''
                plt.savefig(os.path.join(central_limit_kimpap_folder, f'plot_bootstrap_grid_{pattern_class}{animal_suffix}_epoch{figappend}.png'), dpi=300)
                plt.savefig(os.path.join(central_limit_kimpap_folder, f'plot_bootstrap_grid_{pattern_class}{animal_suffix}_epoch{figappend}.pdf'))
                

    # The function is ready for use once the data is reloaded.
    base_query_str = 'pattern_cca1 == 2 & pattern_cca2==10'
    generate_and_save_epoch_bootstrap_plots(df_clean, base_query_str, ['theta', 'delta', 'ripple'], animals=['all'], combine_epochs=True)
    pattern_class_names = ['theta', 'delta', 'ripple']
    plt.close('all')

    generate_and_save_epoch_bootstrap_plots(df_clean, base_query_str, ['theta', 'delta', 'ripple'], animals=list(df_clean.animal.unique()), combine_epochs=True)
    plt.close('all')

    base_query_str = 'pattern_cca1 == 2 & pattern_cca2==0'
    generate_and_save_epoch_bootstrap_plots(df_clean, base_query_str, ['theta', 'delta', 'ripple'], animals=['all'], combine_epochs=True, figappend="_match")
    plt.close('all')

    generate_and_save_epoch_bootstrap_plots(df_clean, base_query_str, ['theta', 'delta', 'ripple'], animals=list(df_clean.animal.unique()), combine_epochs=True, figappend="_match")
    plt.close('all')

