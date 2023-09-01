# Assumes eventuv_plots.py

# PLOT: Bar plot of number of unique animals per uv_component, animal, and genH 

# Group by 'uv_component', 'animal', and 'genH' and count unique animals
def plot_animal_counts(df):
    grouped = df.groupby(['uv_components', 'genH']).animal.nunique().reset_index()
    grouped.rename(columns={'animal': 'animal_count'}, inplace=True)

    # Create the bar plot
    plt.figure(figsize=(15, 6))
    sns.barplot(x='uv_components', y='animal_count', hue='genH', data=grouped)

    plt.title('Number of Unique Animals per uv_component, animal, and genH')
    plt.ylabel('Number of Unique Animals')
    plt.xlabel('UV Component')
    plt.legend(title='genH')
    plt.show()

plot_animal_counts(df_clean)
