# Assumes eventuv_plots.py

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import mplcyberpunk

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

# PLOT: Let's figure out the distribution of 

g=sns.pairplot(data=df_clean, vars=['pattern_cca1', 'pattern_cca2', 'patterns', 'genH', 'highlow','epoch'])
# add axis spines on all sides
for ax in g.axes.flatten():
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    # set colors to black
    ax.spines['top'].set_color('black')
    ax.spines['right'].set_color('black')
    # set linewidths to 2
    ax.spines['top'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
# mdf_clean.pattern_cca1.unique()
# g.map_lower(lambda *pos,**kws : plt.cla())
mplcyberpunk.make_lines_glow()
plt.pause(0.1)
plt.savefig(os.path.join(figfolder, 'plot_pairplot_patterns_and_genH.png'))
plt.savefig(os.path.join(figfolder, 'plot_pairplot_patterns_and_genH.pdf'))


df_matrix.reset_index(inplace=True)

g=sns.pairplot(data=df_matrix, vars=['pattern_cca1', 'pattern_cca2', 'patterns', 'genH', 'highlow', 'epoch'])
# mdf_clean.pattern_cca1.unique()
mplcyberpunk.make_lines_glow()
g.map_lower(lambda *pos,**kws : plt.cla())
plt.suptitle("Pairplot of patterns, genH, and highlow - df_matrix")
plt.pause(0.1)

g=sns.pairplot(data=df_matrix, vars=['pattern_cca1', 'pattern_cca2', 'patterns', 'genH', 'highlow'])
# mdf_clean.pattern_cca1.unique()
mplcyberpunk.make_lines_glow()
g.map_lower(lambda *pos,**kws : plt.cla())
plt.suptitle("Pairplot of patterns, genH, and highlow - df_matrix")
plt.pause(0.1)


df_clean.groupby('patterns').pattern_cca2.value_counts().sort_index()

df_clean.groupby(['genH', 'highlow']).animal.value_counts()

df_clean.groupby('genH_highlow').highlow.value_counts()
