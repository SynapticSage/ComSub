from scipy import stats
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
plt.close('all')


#--------------------

# For theta, hpc-pfc, let's grab the u and v scatter plot measured against the overall pattern
plt.close('all')
def qfilt(df, X):
    lower_quantile = df[X].quantile(0.01)
    upper_quantile = df[X].quantile(0.99)
    filtered_df = df[(df[X] >= lower_quantile) & (df[X] <= upper_quantile)]
    return filtered_df
tmp=df_clean.query('highlow == "high" & pattern_cca1 == 2 & pattern_cca2==10 & uv_components == 1').groupby('genH').sample(1_000)
tmp = qfilt(qfilt(tmp, 'event_u_values'), 'event_v_values')
g=sns.relplot(data=tmp, kind='scatter', 
                 x='event_u_values', y='event_v_values', alpha=1, legend=False,
                hue="pattern_class_name", col="pattern_class_name", row="genH")
g.map(lambda *pos,**kws: plt.axvline(0, color='black', linestyle='--'))
g.map(lambda *pos,**kws: plt.axhline(0, color='black', linestyle='--'))
plt.show()
plt.savefig(os.path.join(figfolder, f'plot_scatter_eventuv_hpcpfc.png'))
plt.savefig(os.path.join(figfolder, f'plot_scatter_eventuv_hpcpfc.pdf'))

# tmp=df_clean.query('highlow == "high" & pattern_cca1 == 1 & pattern_cca2==10').groupby('genH').sample(1_000)
# tmp = qfilt(qfilt(tmp, 'event_u_values'), 'event_v_values')
# g=sns.relplot(data=tmp, kind='scatter', 
#                  x='event_u_values', y='event_v_values', alpha=1, legend=False,
#                 hue="pattern_class_name", col="pattern_class_name", row="genH")
# g.map(lambda *pos,**kws: plt.axvline(0, color='black', linestyle='--'))
# g.map(lambda *pos,**kws: plt.axhline(0, color='black', linestyle='--'))
# plt.show()
# plt.savefig(os.path.join(figfolder, f'plot_scatter_eventuv_hpchpc.png'))
# plt.savefig(os.path.join(figfolder, f'plot_scatter_eventuv_hpchpc.pdf'))

tmp=df_clean.query('highlow == "high" & pattern_cca1 == 2 & pattern_cca2==10 & uv_components==1').groupby('genH').sample(1_000)
tmp = qfilt(qfilt(tmp, 'projection_score'), 'perpendicular_score')
g=sns.relplot(data=tmp, kind='scatter', 
                 x='projection_score', y='perpendicular_score', alpha=1, legend=False,
                hue="pattern_class_name", col="pattern_class_name", row="genH")
g.map(lambda *pos,**kws: plt.axvline(0, color='black', linestyle='--'))
g.map(lambda *pos,**kws: plt.axhline(0, color='black', linestyle='--'))
for ax in g.axes.flat:
    # find max abs x and max abs y
    xmax = max(abs(ax.get_xlim()[0]), abs(ax.get_xlim()[1]))
    ymax = max(abs(ax.get_ylim()[0]), abs(ax.get_ylim()[1]))
    mx = max(xmax, ymax)
    ax.set_xlim(-mx, mx)
    ax.set_ylim(-mx, mx)
plt.show()
plt.savefig(os.path.join(figfolder, f'plot_scatter_projperp_hpcpfc.png'))
plt.savefig(os.path.join(figfolder, f'plot_scatter_projperp_hpcpfc.pdf'))


# Let's do the above plot, but a kdeplot, sample=10_000
tmp=df_clean.query('pattern_cca1 == 2 & pattern_cca2==10 & uv_components == 1').groupby('genH').sample(7_500)
tmp = qfilt(qfilt(tmp, 'projection_score'), 'perpendicular_score')
g=sns.FacetGrid(data=tmp, 
                row="pattern_class_name", col="genH_highlow", hue="pattern_class_name")
g.map(sns.kdeplot, "projection_score", "perpendicular_score", alpha=0.5, shade=True)
g.map(sns.scatterplot, "projection_score", "perpendicular_score", alpha=0.02)
g.map(lambda *pos,**kws: plt.axvline(0, color='black', linestyle='--'))
g.map(lambda *pos,**kws: plt.axhline(0, color='black', linestyle='--'))
for ax in g.axes.flat:
    # find max abs x and max abs y
    xmax = max(abs(ax.get_xlim()[0]), abs(ax.get_xlim()[1]))
    ymax = max(abs(ax.get_ylim()[0]), abs(ax.get_ylim()[1]))
    mx = max(xmax, ymax)
    ax.set_xlim(-mx, mx)
    ax.set_ylim(-mx, mx)
    # Exchange | with \n in title
    title = ax.get_title().replace('|', '\n')
    ax.set_title(title)
plt.savefig(os.path.join(figfolder, f'plot_kde_projperp_hpcpfc.png'))
plt.savefig(os.path.join(figfolder, f'plot_kde_projperp_hpcpfc.pdf'))

# A projection and perpendicular axis version of this
for i in range(1,4):
    # Let's do the above plot, but a kdeplot, sample=20_000
    tmp=(df_clean.query('pattern_cca1 == 2 & pattern_cca2==10 & (highlow == "high" | highlow == "low") & uv_components == @i')
         .groupby('genH').sample(10_000))
    tmp = qfilt(qfilt(tmp, 'projection_score'), 'perpendicular_score')
    g=sns.FacetGrid(data=tmp, 
                    row="pattern_class_name", col="genH", hue="highlow")
    g.map(sns.kdeplot, "projection_score", "perpendicular_score", alpha=0.5, shade=True)
    g.map(sns.scatterplot, "projection_score", "perpendicular_score", alpha=0.02)
    g.map(lambda *pos,**kws: plt.axvline(0, color='black', linestyle='--'))
    g.map(lambda *pos,**kws: plt.axhline(0, color='black', linestyle='--'))
    for ax in g.axes.flat:
        # find max abs x and max abs y
        xmax = max(abs(ax.get_xlim()[0]), abs(ax.get_xlim()[1]))
        ymax = max(abs(ax.get_ylim()[0]), abs(ax.get_ylim()[1]))
        mx = max(xmax, ymax)
        ax.set_xlim(-mx, mx)
        ax.set_ylim(-mx, mx)
        # Exchange | with \n in title
        title = ax.get_title().replace('|', '\n')
        ax.set_title(title)
    # remove alpha from the points for the legend
    g.add_legend()
    # Adjust alpha and size of legend markers
    if g._legend is not None:
        for lh in g._legend.legendHandles:
            lh.set_alpha(1)
            lh._sizes = [50]  # this attribute is a list
    plt.suptitle(f'UV Component {i} - HPC-PFC')
    plt.subplots_adjust(top=0.9)
    plt.savefig(os.path.join(figfolder, f'plot_kde_projperp_hpcpfc_hue=highlow_uvcomp={i}.png'))
    plt.savefig(os.path.join(figfolder, f'plot_kde_projperp_hpcpfc_hue=highlow_uvcomp={i}.pdf'))

# An event U and V version of this
for i in range(1,4):
    # Let's do the above plot, but a kdeplot, sample=20_000
    tmp=(df_clean.query('pattern_cca1 == 2 & pattern_cca2==10 & (highlow == "high" | highlow == "low") & uv_components == @i')
         .groupby('genH').sample(10_000))
    tmp = qfilt(qfilt(tmp, 'event_u_values'), 'event_v_values')
    g=sns.FacetGrid(data=tmp, 
                    row="pattern_class_name", col="genH", hue="highlow")
    g.map(sns.kdeplot, "event_u_values", "event_v_values", alpha=0.5, shade=True)
    g.map(sns.scatterplot, "event_u_values", "event_v_values", alpha=0.02)
    g.map(lambda *pos,**kws: plt.axvline(0, color='black', linestyle='--'))
    g.map(lambda *pos,**kws: plt.axhline(0, color='black', linestyle='--'))
    for ax in g.axes.flat:
        # find max abs x and max abs y
        xmax = max(abs(ax.get_xlim()[0]), abs(ax.get_xlim()[1]))
        ymax = max(abs(ax.get_ylim()[0]), abs(ax.get_ylim()[1]))
        mx = max(xmax, ymax)
        ax.set_xlim(-mx, mx)
        ax.set_ylim(-mx, mx)
        # Exchange | with \n in title
        title = ax.get_title().replace('|', '\n')
        ax.set_title(title)
    # remove alpha from the points for the legend
    g.add_legend()
    # Adjust alpha and size of legend markers
    if g._legend is not None:
        for lh in g._legend.legendHandles:
            lh.set_alpha(1)
            lh._sizes = [50]  # this attribute is a list
    plt.suptitle(f'UV Component {i} - HPC-PFC')
    plt.subplots_adjust(top=0.9)
    plt.savefig(os.path.join(figfolder, f'plot_kde_eventuv_hpcpfc_hue=highlow_uvcomp={i}.png'))
    plt.savefig(os.path.join(figfolder, f'plot_kde_eventuv_hpcpfc_hue=highlow_uvcomp={i}.pdf'))


# Let's do a jointplot of the projection and perpendicular scores for theta, hpc-pfc, uvcomp=1
for i in range(1,4):
    tmp=(df_clean.query('pattern_cca1 == 2 & pattern_cca2==10 & uv_components == 1 & (genH_highlow == "coherence_high" | genH_highlow == "power_high")')
    .groupby('genH_highlow').sample(5_000))
    assert(tmp.genH_highlow.unique().size == 2)
    tmp = qfilt(qfilt(tmp, 'projection_score'), 'perpendicular_score')
    g=sns.JointGrid(data=tmp, x="projection_score", y="perpendicular_score", 
                    hue="genH_highlow")
    g.plot_marginals(sns.histplot, alpha=0.5, multiple="stack")
    g.plot_joint(sns.kdeplot, alpha=0.5, shade=True)
    g.ax_joint.set_xlim(-5, 5)
    g.ax_joint.set_ylim(-5, 5)
    # Test for equal projection distribution (ks-test)
    ks, p = stats.ks_2samp(tmp.query('genH_highlow == "coherence_high"').projection_score,
                           tmp.query('genH_highlow == "power_high"').projection_score)
    ks2, p2 = stats.ks_2samp(tmp.query('genH_highlow == "coherence_high"').perpendicular_score,
                            tmp.query('genH_highlow == "power_high"').perpendicular_score)
    # Two sample t-tests
    t, pt = stats.ttest_ind(tmp.query('genH_highlow == "coherence_high"').projection_score,
                            tmp.query('genH_highlow == "power_high"').projection_score)
    t2, pt2 = stats.ttest_ind(tmp.query('genH_highlow == "coherence_high"').perpendicular_score,
                            tmp.query('genH_highlow == "power_high"').perpendicular_score)
    # KS-test values, and for ttest, let's show the difference and the pvalue
    g.fig.suptitle(f'UV Component {i} - HPC-PFC - KS test proj={p:.3f} perp={p2:.3f}\n'
                   f'T-test proj={t:.3f} ({pt:.3f}) perp={t2:.3f} ({pt2:.3f})')
    g.fig.subplots_adjust(top=0.9)
    # let's create a very descriptive title
    plt.savefig(os.path.join(f"{figfolder}", f'plot_joint_projperp_hpcpfc_hue=genH_highlow_uvcomp={i}.png'))
    plt.savefig(os.path.join(f"{figfolder}", f'plot_joint_projperp_hpcpfc_hue=genH_highlow_uvcomp={i}.pdf'))


# Let's do high coherence versus low-coherence, same
for i in range(1,4):
    tmp=(df_clean.query('pattern_cca1 == 2 & pattern_cca2==10 & uv_components == 1 & (genH_highlow == "coherence_high" | genH_highlow == "coherence_low")')
    .groupby('genH_highlow').sample(5_000))
    assert(tmp.genH_highlow.unique().size == 2)
    tmp = qfilt(qfilt(tmp, 'projection_score'), 'perpendicular_score')
    g=sns.JointGrid(data=tmp, x="projection_score", y="perpendicular_score", 
                    hue="genH_highlow")
    g.plot_marginals(sns.histplot, alpha=0.5, multiple="stack")
    g.plot_joint(sns.kdeplot, alpha=0.5, shade=True)
    g.ax_joint.set_xlim(-5, 5)
    g.ax_joint.set_ylim(-5, 5)
    # Test for equal projection distribution (ks-test)
    ks, p = stats.ks_2samp(tmp.query('genH_highlow == "coherence_high"').projection_score,
                           tmp.query('genH_highlow == "coherence_low"').projection_score)
    ks2, p2 = stats.ks_2samp(tmp.query('genH_highlow == "coherence_high"').perpendicular_score,
                            tmp.query('genH_highlow == "coherence_low"').perpendicular_score)
    # Two sample t-tests
    t, pt = stats.ttest_ind(tmp.query('genH_highlow == "coherence_high"').projection_score,
                            tmp.query('genH_highlow == "coherence_low"').projection_score)
    t2, pt2 = stats.ttest_ind(tmp.query('genH_highlow == "coherence_high"').perpendicular_score,
                            tmp.query('genH_highlow == "coherence_low"').perpendicular_score)
    # KS-test values, and for ttest, let's show the difference and the pvalue
    g.fig.suptitle(f'UV Component {i} - HPC-PFC - KS test proj={p:.3f} perp={p2:.3f}\n'
                   f'T-test proj={t:.3f} ({pt:.3f}) perp={t2:.3f} ({pt2:.3f})')
    g.fig.subplots_adjust(top=0.9)
    # let's create a very descriptive title
    plt.savefig(os.path.join(f"{figfolder}", f'coherence_plot_joint_projperp_hpcpfc_hue=genH_highlow_uvcomp={i}.png'))
    plt.savefig(os.path.join(f"{figfolder}", f'coherence_plot_joint_projperp_hpcpfc_hue=genH_highlow_uvcomp={i}.pdf'))


# Let's do the same for power high and power low now
plt.close('all')
for pattern_class_name in ['theta', 'delta', 'ripple']:
    for i in range(1,4):
        tmp=(df_clean.query('pattern_cca1 == 2 & pattern_cca2==10 & uv_components == @i & (genH_highlow == "power_high" | genH_highlow == "power_low") & pattern_class_name == @pattern_class_name')
        .groupby('genH_highlow').sample(2_500))
        assert(tmp.genH_highlow.unique().size == 2)
        tmp = qfilt(qfilt(tmp, 'projection_score'), 'perpendicular_score')
        g=sns.JointGrid(data=tmp, x="projection_score", y="perpendicular_score", 
                        hue="genH_highlow")
        g.plot_marginals(sns.histplot, alpha=0.5, multiple="stack")
        g.plot_joint(sns.kdeplot, alpha=0.5, shade=True)
        mx = tmp[['projection_score', 'perpendicular_score']].abs().max().max()
        g.ax_joint.set_xlim(-mx, mx)
        g.ax_joint.set_ylim(-mx, mx)
        # Test for equal projection distribution (ks-test)
        ks, p = stats.ks_2samp(tmp.query('genH_highlow == "power_high"').projection_score,
                               tmp.query('genH_highlow == "power_low"').projection_score)
        ks2, p2 = stats.ks_2samp(tmp.query('genH_highlow == "power_high"').perpendicular_score,
                                tmp.query('genH_highlow == "power_low"').perpendicular_score)
        # Two sample t-tests
        t, pt = stats.ttest_ind(tmp.query('genH_highlow == "power_high"').projection_score,
                                tmp.query('genH_highlow == "power_low"').projection_score)
        t2, pt2 = stats.ttest_ind(tmp.query('genH_highlow == "power_high"').perpendicular_score,
                                tmp.query('genH_highlow == "power_low"').perpendicular_score)
        # KS-test values, and for ttest, let's show the difference and the pvalue
        g.fig.suptitle(f'{pattern_class_name} - UV Component {i} - HPC-PFC - KS test proj={p:.3f} perp={p2:.3f}\n' 
                f'UV Component {i} - HPC-PFC - KS test proj={p:.3f} perp={p2:.3f}\n' +
                f'T-test proj={t:.3f} ({pt:.3f}) perp={t2:.3f} ({pt2:.3f})')
        g.fig.subplots_adjust(top=0.9)
        # let's create a very descriptive title
        plt.savefig(os.path.join(f"{figfolder}", f'power_plot_joint_projperp_hpcpfc_hue=genH_highlow_uvcomp={i}_{pattern_class_name}.png'))
        plt.savefig(os.path.join(f"{figfolder}", f'power_plot_joint_projperp_hpcpfc_hue=genH_highlow_uvcomp={i}_{pattern_class_name}.pdf'))

# Let's do the above plot, but a kdeplot, sample=10_000
for pattern_class_name in ['theta', 'delta', 'ripple']:
    tmp=df_clean.query('pattern_cca1 == 2 & pattern_cca2==10 & pattern_class_name == @pattern_class_name').groupby('genH').sample(7_500)
    tmp = qfilt(qfilt(tmp, 'projection_score'), 'perpendicular_score')
    colors_half1 = cm.Oranges(np.linspace(0, 1, 3))
    colors_half2 = cm.Blues(np.linspace(0, 1, 3))
    pal = sns.color_palette([*colors_half1, *colors_half2])
    g=sns.FacetGrid(data=tmp, row="animal", hue="genH_highlow", col="genH_highlow",
                    col_order=['power_low', 'power_mid', 'power_high', 'coherence_low', 'coherence_mid', 'coherence_high'], 
                    hue_order=['power_low', 'power_mid', 'power_high', 'coherence_low', 'coherence_mid', 'coherence_high'],
                    palette=pal)
    g.map(sns.kdeplot, "projection_score", "perpendicular_score", alpha=0.5, shade=True)
    g.map(sns.scatterplot, "projection_score", "perpendicular_score", alpha=0.02)
    g.map(lambda *pos,**kws: plt.axvline(0, color='black', linestyle='--'))
    g.map(lambda *pos,**kws: plt.axhline(0, color='black', linestyle='--'))
    for ax in g.axes.flat:
        # find max abs x and max abs y
        xmax = max(abs(ax.get_xlim()[0]), abs(ax.get_xlim()[1]))
        ymax = max(abs(ax.get_ylim()[0]), abs(ax.get_ylim()[1]))
        mx = max(xmax, ymax)
        ax.set_xlim(-mx, mx)
        ax.set_ylim(-mx, mx)
        # Exchange | with \n in title
        title = ax.get_title().replace('|', '\n')
        ax.set_title(title)
    plt.savefig(os.path.join(figfolder, f'{pattern_class_name}_animal_plot_kde_projperp_hpcpfc.png'))
    plt.savefig(os.path.join(figfolder, f'{pattern_class_name}_animal_plot_kde_projperp_hpcpfc.pdf'))


# Let's do the above plot, but a kdeplot, sample=10_000
for pattern_class_name in ['theta', 'delta', 'ripple']:
    for comp in range(1,6):
        tmp=df_clean.query('pattern_cca1 == 2 & pattern_cca2==10 & pattern_class_name == @pattern_class_name & uv_components == @comp').groupby('genH').sample(7_500)
        tmp = qfilt(qfilt(tmp, 'projection_score'), 'perpendicular_score')
        colors_half1 = cm.Oranges(np.linspace(0, 1, 3))
        colors_half2 = cm.Blues(np.linspace(0, 1, 3))
        pal = sns.color_palette([*colors_half1, *colors_half2])
        g=sns.FacetGrid(data=tmp, row="animal", hue="genH_highlow", col="genH_highlow",
                        col_order=['power_low', 'power_mid', 'power_high', 'coherence_low', 'coherence_mid', 'coherence_high'], 
                        hue_order=['power_low', 'power_mid', 'power_high', 'coherence_low', 'coherence_mid', 'coherence_high'],
                        palette=pal)
        g.map(sns.kdeplot, "projection_score", "perpendicular_score", alpha=0.5, shade=True)
        g.map(sns.scatterplot, "projection_score", "perpendicular_score", alpha=0.02)
        g.map(lambda *pos,**kws: plt.axvline(0, color='black', linestyle='--'))
        g.map(lambda *pos,**kws: plt.axhline(0, color='black', linestyle='--'))
        for ax in g.axes.flat:
            # find max abs x and max abs y
            xmax = max(abs(ax.get_xlim()[0]), abs(ax.get_xlim()[1]))
            ymax = max(abs(ax.get_ylim()[0]), abs(ax.get_ylim()[1]))
            mx = max(xmax, ymax)
            ax.set_xlim(-mx, mx)
            ax.set_ylim(-mx, mx)
            # Exchange | with \n in title
            title = ax.get_title().replace('|', '\n')
            ax.set_title(title)
        plt.savefig(os.path.join(figfolder, f'{pattern_class_name}_animal_plot_kde_projperp_hpcpfc_comp={comp}.png'))
        plt.savefig(os.path.join(figfolder, f'{pattern_class_name}_animal_plot_kde_projperp_hpcpfc_comp={comp}.pdf'))

