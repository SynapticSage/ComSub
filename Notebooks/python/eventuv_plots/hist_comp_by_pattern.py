# run eventuv_plots.py first

# PLOT: Hist plot of on/off commsub -------------
# Create a subplot grid with one row for each uv_component and one column for each pattern
def plot_hist(df_clean, thing='on_commsub_mag', pattern_cca2=overall_pattern, **kws):
    fig, axs = plt.subplots(len(uv_components), len(patterns), figsize=(4*len(patterns), 4*len(uv_components)), sharey=True)
    if len(uv_components) == 1:
        axs = [axs]
    if len(patterns) == 1:
        axs = [[ax] for ax in axs]
    
    df_clean = df_clean.sort_values(['patterns'])
    q = df_clean[thing].quantile(0.99)
    qlow = df_clean[thing].quantile(0.001)
    
    for i, uv_component in enumerate(uv_components):
        for j, pattern in enumerate(patterns):
            ax = axs[i][j]
            
            # Check if pattern_cca2 is set to 'match'
            if pattern_cca2 == 'match':
                pattern_cca2_value = pattern
            else:
                pattern_cca2_value = pattern_cca2
                
            for k, genH in enumerate(genH_highlow):
                df_sub = df_clean[
                    (df_clean['uv_components'] == uv_component) & 
                    (df_clean['patterns'] == pattern) & 
                    (df_clean['genH_highlow'] == genH) & 
                    (df_clean['pattern_cca2'] == pattern_cca2_value)
                ]
                if all(df_sub.empty) or all(df_sub[thing].isna()):
                    print(f'No data for uv_component {uv_component}, pattern {pattern}, genH {genH}, pattern_cca2 {pattern_cca2_value}')
                    continue
                
                q_75th = df_sub[thing].quantile(0.75)
                q_25th = df_sub[thing].quantile(0.25)
                q_90 = df_sub[thing].mean()
                
                g = sns.histplot(df_sub[thing], 
                                 ax=ax, kde=False, label=genH, binrange=(-10, 10), **kws)
                palette = sns.color_palette()
                color = palette[k % len(palette)]
                ax.axvline(q_25th, color=color, linestyle='--', label=f'75th percentile {genH}')
                ax.axvline(q_75th, color=color, linestyle='--', label=f'75th percentile {genH}')
                ax.axvline(q_90, color=color, linestyle=':')
                
            if i == 0:
                ax.set_title(f'Pattern {pattern}')
            if j == 0:
                ax.set_ylabel(f'UV Component {int(uv_component)}')
            
            if i == 0 and (pattern == 0 or pattern == 3):
                ax.legend(title='genH')
            ax.set_xlim(qlow, q)
            
    plt.tight_layout()
    plt.suptitle(f'{thing}')
    plt.subplots_adjust(top=0.95)
    plt.show()
    return g



plot_hist(df_clean, thing='on_commsub', stat="density", alpha=0.3)
plt.savefig(os.path.join(figfolder, 'plot_hist_thing=on_commsub.png'))
plt.savefig(os.path.join(figfolder, 'plot_hist_thing=on_commsub.pdf'))

plot_hist(df_clean, thing='on_commsub_mag', stat="density", alpha=0.3)
plt.savefig(os.path.join(figfolder, 'plot_hist_thing=on_commsub_mag.png'))
plt.savefig(os.path.join(figfolder, 'plot_hist_thing=on_commsub_mag.pdf'))
g=plot_hist(df_clean, thing='on_commsub_mag', stat="density", cumulative=True, alpha=0.45)
plt.savefig(os.path.join(figfolder, 'plot_hist_thing=on_commsub_mag_cumulative.png'))
plt.savefig(os.path.join(figfolder, 'plot_hist_thing=on_commsub_mag_cumulative.pdf'))
plot_hist(df_clean, thing='scaled_on_commsub_mag')
plt.savefig(os.path.join(figfolder, 'plot_hist_thing=scaled_on_commsub_mag.png'))
plt.savefig(os.path.join(figfolder, 'plot_hist_thing=scaled_on_commsub_mag.pdf'))

plot_hist(df_clean, thing='on_commsub_mag', pattern_cca2='match', stat="density", alpha=0.3)
plt.savefig(os.path.join(figfolder, 'plot_hist_thing=on_commsub_mag_pattern_cca2=match.png'))
plt.savefig(os.path.join(figfolder, 'plot_hist_thing=on_commsub_mag_pattern_cca2=match.pdf'))
plot_hist(df_clean, thing='on_commsub_mag', pattern_cca2='match', stat="density", cumulative=True, alpha=0.45)
plt.savefig(os.path.join(figfolder, 'plot_hist_thing=on_commsub_mag_pattern_cca2=match_cumulative.png'))
plt.savefig(os.path.join(figfolder, 'plot_hist_thing=on_commsub_mag_pattern_cca2=match_cumulative.pdf'))

plot_hist(df_clean, thing='on_over_off', stat="density", alpha=0.3)

os.system('pushover-cli "Finished out plot_hist"')

# PLOT: Hist plot of event_u_values vs event_v_values -------------

def plot_hist_aggregated(df, thing='on_commsub_mag', pattern_cca2=7):
    """
    Create a grid of histograms for each combination of pattern, colored by 'genH'.
    Each histogram represents the distribution of the 'thing' column values.
    A vertical dashed line is drawn at the 75th percentile of the 'thing' column for each 'genH'.
    """
    # Identify the unique patterns and genH values
    patterns = sorted(df['patterns'].unique())
    genH_values = df['genH'].unique()

    # Create a subplot grid with one column for each pattern
    fig, axs = plt.subplots(1, len(patterns), figsize=(5 * len(patterns), 5), sharey=True)
    # Adjust for case when there is only one pattern
    if len(patterns) == 1:
        axs = [axs]
    # Calculate the 99th percentile of the 'thing' column for the x-axis limit
    q = df[thing].quantile(0.99)
    # For each pattern
    for j, pattern in enumerate(sorted(patterns)):
        ax = axs[j]
        # For each genH value
        for k, genH in enumerate(genH_highlow):
            # Filter the DataFrame for the current pattern, genH, and pattern_cca2
            df_sub = df[(df['patterns'] == pattern) & (df['genH_highlow'] == genH) & (df['pattern_cca2'] == pattern_cca2)]
            if df_sub.empty:
                print(f'No data for pattern {pattern}, genH {genH}, pattern_cca2 {pattern_cca2}')
                continue
            # Calculate the 75th percentile of the 'thing' column
            q_half = df_sub[thing].quantile(0.75)
            # Plot a histogram of the 'thing' column values
            g = sns.histplot(df_sub[thing], ax=ax, kde=False, label=genH,
                             binrange=(0, 10), stat='density',
                             common_norm=False)
            # Draw a vertical dashed line at the 75th percentile
            palette = sns.color_palette()
            color = palette[k % len(palette)]
            ax.axvline(q_half, color=color, linestyle='--', label=f'75th percentile {genH}')
        # Set titles and labels
        ax.set_title(f'Pattern {pattern}')
        if j == 0:
            ax.set_ylabel(f'UV Component Aggregated')
        ax.legend(title='genH')
        ax.set_xlim(0, q)

    plt.tight_layout()
    plt.suptitle(f'{thing}')
    plt.subplots_adjust(top=0.95)
    plt.show()

# Test the function with the current DataFrame
plot_hist_aggregated(df_clean, thing='on_commsub_mag')
plt.savefig(os.path.join(figfolder, 'plot_hist_aggregated_thing=on_commsub_mag.png'))
plot_hist_aggregated(df_clean, thing='scaled_on_commsub_mag')
plt.savefig(os.path.join(figfolder, 'plot_hist_aggregated_thing=scaled_on_commsub_mag.png'))

