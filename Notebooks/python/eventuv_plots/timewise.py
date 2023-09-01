
# ========
# ALERT:
# UV magnitude OVER TIME 󰥔
# ========

# TODO: CUT OUT DIFF GREATER THAN 1 minute

df_matrix = prep_uv_melt(df_clean)
df_matrix = prep_uv_magnitude_over_time(df_matrix)
df_matrix.query('pattern_cca2 == @overall_pattern', inplace=True)
df_matrix = df_matrix.reset_index()

#9 Melt the data for easy plotting
df_melted = df_matrix.melt(id_vars=['epoch','time'], value_vars=['magnitude_u', 'magnitude_v'], var_name='component', value_name='magnitude')
# Plot using relplot
print("Plotting...")
g = sns.relplot(data=df_melted, x='time', y='magnitude', kind='line', hue='component', height=5, aspect=2, col='epoch', col_wrap=4)
g.set_axis_labels("Event Time", "Magnitude")
g.tight_layout()
plt.show()
plt.close()


df_melted = df_matrix.melt(id_vars=['time','epoch','animal','genH'], value_vars=['magnitude_u_smooth', 'magnitude_v_smooth'], var_name='component', value_name='magnitude')
# Plot using relplot
print("Plotting...")
g = sns.relplot(data=df_melted, x='time', y='magnitude', kind='line',
                hue='component', height=5, aspect=2, alpha=0.5, col='epoch',
                row='animal', facet_kws=dict(sharex=False))
g.set_axis_labels("Event Time", "Magnitude")
g.tight_layout()
plt.show()
plt.savefig(os.path.join(figfolder,'match_to_genH_overtime_hue=component_row=animal_col=epoch.png'), dpi=300)
plt.savefig(os.path.join(figfolder,'match_to_genH_overtime_hue=component_row=animal_col=epoch.pdf'), dpi=300)

g = sns.relplot(data=df_melted, x='time', y='magnitude', kind='line',
                hue='genH', height=5, aspect=2, alpha=0.5, col='epoch',
                row='animal', facet_kws=dict(sharex=False))
g.set_axis_labels("Event Time", "Magnitude")
g.tight_layout()
plt.show()
plt.savefig(os.path.join(figfolder,'match_to_genH_overtime_hue=component_row=animal_col=epoch.png'), dpi=300)
plt.savefig(os.path.join(figfolder,'match_to_genH_overtime_hue=component_row=animal_col=epoch.pdf'), dpi=300)


# Plot using relplot
print("Plotting...")
g = sns.relplot(data=df_matrix, x='time_bin', y='magnitude_u', kind='line', hue='genH', height=5, aspect=2)
g.set_axis_labels("Time Bin", "Magnitude")
g.tight_layout()
plt.show()


# Plot using relplot
print("Plotting...")
g = sns.relplot(data=df_matrix, x='time_bin', y='magnitude_u', kind='line', hue='genH_highlow', height=5, aspect=2)
g.set_axis_labels("Time Bin", "Magnitude")
g.tight_layout()
plt.show()



# Plot using relplot
print("Plotting...")
g = sns.relplot(data=df_matrix, x='time_bin', y='magnitude_v', kind='line', hue='genH', height=5, aspect=2)
g.set_axis_labels("Time Bin", "Magnitude")
g.tight_layout()
plt.show()


# Plot using relplot
print("Plotting...")
g = sns.relplot(data=df_matrix, x='time_bin', y='magnitude_v', kind='line', hue='genH_highlow', height=5, aspect=2)
g.set_axis_labels("Time Bin", "Magnitude")
g.tight_layout()
plt.show()


# Plot using relplot
print("Plotting...")
g = sns.relplot(data=df_matrix.query('genH == "coherence"'), x='time_bin', y='magnitude_u', kind='line', hue='genH_highlow', height=5, aspect=2)
g.set_axis_labels("Time Bin", "Magnitude")
g.tight_layout()
plt.show()
