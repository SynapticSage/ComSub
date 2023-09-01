# Description: This script is used to generate scatter plots of the event u and v values
# Runs in the eventuv_plots.py script


# PLOT: : IMPORTANT
import seaborn as sns
import matplotlib.pyplot as plt
# Create a FacetGrid that arranges the data by uv_components (rows) and patterns (columns)
g = sns.FacetGrid(df_clean, col='patterns', row='uv_components', hue='genH',
                  height=4, aspect=1, palette='viridis')
# Apply a scatterplot onto each facet with increased alpha transparency
g.map(sns.scatterplot, 'event_u_values', 'event_v_values', alpha=0.05)
# Add a diagonal line, text, and vertical/horizontal lines at x=0 and y=0 for each subplot
def add_lines(x, y, **kwargs):
    ax = plt.gca()
    ax.axline((0, 0), slope=1, color='black', linestyle='--')
    ax.text(0.6, 0.7, 'comm sub', transform=ax.transAxes, color='black')
    ax.axvline(0, color='black', linestyle=':', linewidth=0.5)
    ax.axvline(5, color='black', linestyle=':', linewidth=0.5)
    ax.axvline(-5, color='black', linestyle=':')
    ax.axhline(-5, color='black', linestyle=':', linewidth=0.5)
    ax.axhline(5, color='black', linestyle=':', linewidth=0.5)
    ax.axhline(0, color='black', linestyle=':')

g.map(add_lines, 'event_u_values', 'event_v_values')
# Add a legend
g.add_legend()
# Display the plot
plt.show()
plt.savefig(os.path.join(figfolder, 'plot_scatter_eventuv.png'))
plt.savefig(os.path.join(figfolder, 'plot_scatter_eventuv.pdf'))

# PLOT: : IMPORTANT
import seaborn as sns
import matplotlib.pyplot as plt
# Create a FacetGrid that arranges the data by uv_components (rows) and patterns (columns)
g = sns.FacetGrid(df_clean, col='highlow', row='uv_components', hue='genH',
                  height=4, aspect=1, palette='viridis')
# Apply a scatterplot onto each facet with increased alpha transparency
g.map(sns.scatterplot, 'event_u_values', 'event_v_values', alpha=0.05)
# Add a diagonal line, text, and vertical/horizontal lines at x=0 and y=0 for each subplot
def add_lines(x, y, **kwargs):
    ax = plt.gca()
    ax.axline((0, 0), slope=1, color='black', linestyle='--')
    ax.text(0.6, 0.7, 'comm sub', transform=ax.transAxes, color='black')
    ax.axvline(0, color='black', linestyle=':', linewidth=0.5)
    ax.axvline(5, color='black', linestyle=':', linewidth=0.5)
    ax.axvline(-5, color='black', linestyle=':')
    ax.axhline(-5, color='black', linestyle=':', linewidth=0.5)
    ax.axhline(5, color='black', linestyle=':', linewidth=0.5)
    ax.axhline(0, color='black', linestyle=':')

g.map(add_lines, 'event_u_values', 'event_v_values')
# Add a legend
g.add_legend()
# Display the plot
plt.show()
plt.savefig(os.path.join(figfolder, 'plot_scatter_eventuv_highlow.png'))
plt.savefig(os.path.join(figfolder, 'plot_scatter_eventuv_highlow.pdf'))

