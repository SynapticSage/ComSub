def plot_bar(data, y_value, query=None):
    plt.figure(figsize=(12, 6))
    if query:
        data = data.query(query)
    g=sns.barplot(data=data, x='patterns', y=y_value, hue='genH_highlow', errorbar=('ci', 99))
    plt.title(f'Bar Plot for {y_value}\n(99% CI)')
    plt.show()

# Function that takes the sqrt(sum of uv_components squared) 
def mag2(X):
    X = np.sqrt(np.sum(X**2, axis=0))
    return X

# Plotting for projection score
plt.close('all')
df_clean_uvgroup = df_clean.groupby(index).agg({'projection_score': mag2,
                                                'perpendicular_score': mag2,
                                                'abs_projection_score': mag2,
                                                'abs_perpendicular_score': mag2,
                                                'proj_over_perp': 'mean',
                                                'abs_proj_over_perp': 'mean'}).reset_index()
plot_bar(df_clean_uvgroup, 'projection_score')
plot_bar(df_clean_uvgroup, 'perpendicular_score')
plot_bar(df_clean, 'abs_projection_score')
plot_bar(df_clean, 'abs_perpendicular_score')
plot_bar(df_clean, 'proj_over_perp')
plot_bar(df_clean, 'abs_proj_over_perp')

def plot_bar_comp(data, y_value):
    sns.catplot(data=data, x='patterns', y=y_value, hue='genH_highlow', errorbar=('ci', 99),
                row='uv_components', col='genH', kind='bar', sharey=True,
                height=2, aspect=1.5)
    plt.suptitle(f'Bar Plot for {y_value}\n(99% CI)')
    plt.show()

# Plotting for projection score
plt.close('all')
plot_bar_comp(df_clean, 'projection_score')
plot_bar_comp(df_clean, 'perpendicular_score')
plot_bar_comp(df_clean, 'abs_projection_score')
plot_bar_comp(df_clean, 'abs_perpendicular_score')

def plot_bar_animal(data, y_value):
    sns.catplot(data=data, x='patterns', y=y_value, hue='highlow_genH',
                errorbar=('ci', 99), row='uv_components', col='animal',
                kind='bar', sharey=False,
                height=1.5, aspect=1.5)
    plt.suptitle(f'Bar Plot for {y_value}\n(99% CI)')
    plt.show()

plt.close('all')
plot_bar_animal(df_clean, 'abs_projection_score')
plot_bar_animal(df_clean, 'abs_perpendicular_score')
plot_bar_animal(df_clean, 'projection_score')
plot_bar_animal(df_clean, 'perpendicular_score')

def plot_bar_anim_summary(data, y_value, comp_max=None):
    if comp_max:
        data = data.query(f'uv_components <= {comp_max}')
    sns.catplot(data=data, x='patterns', y=y_value, hue='genH_highlow',
                errorbar=('ci', 99), col='animal',
                kind='bar', sharey=False,
                height=1.5, aspect=1.5)
    plt.suptitle(f'Bar Plot for {y_value}\n(99% CI)' + \
                 f' for {comp_max} components' if comp_max else '')
    plt.show()

plt.close('all')
plot_bar_anim_summary(df_clean, 'abs_projection_score',    comp_max=5)
plot_bar_anim_summary(df_clean, 'abs_perpendicular_score', comp_max=5)
plot_bar_anim_summary(df_clean, 'projection_score',        comp_max=5)
plot_bar_anim_summary(df_clean, 'perpendicular_score',     comp_max=5)

