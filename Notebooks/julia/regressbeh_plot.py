import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import mplcyberpunk
plt.style.use("cyberpunk")
plt.ion()
# Plotting styles
sns.set_context("paper", font_scale=1.5, rc={"lines.linewidth": 2.5})

# Folder locations
folder = "/Volumes/MATLAB-Drive/Shared/figures/midpattern=true/data/"
plotfolder = os.path.dirname(folder)
plotfolder = os.path.dirname(plotfolder)
plotfolder = os.path.join(plotfolder, "julia_regressbeh")
if not os.path.exists(plotfolder):
    os.makedirs(plotfolder)

# Read the CSV file
Cdf = pd.read_csv(os.path.join(folder, "coefs_df.csv"))

# Filter out rows where Coefficient is "Intercept"
Cdf.query('Coefficient != "Intercept" & (orthostate == "Orthogonal" | orthostate == "Aligned")', inplace=True)

# ------------------------------------

# First plot with y = "Value"
g = sns.catplot(x="Coefficient", y="Value", hue="Animal", col="Behavior", data=Cdf, kind="box", height=4, aspect=2, legend_out=False, col_wrap=3)
g.set_xticklabels(rotation=45)

for ax in g.axes.flat:
    ax.set_ylim(-1, 1)
    ax.axhline(0, color="black", linestyle="--")
    ax.axvline(11.5, color="black", linestyle="--")
    for bar in ax.patches:
        bar.set_edgecolor("none")
plt.savefig(os.path.join(plotfolder, "coefficients_boxplots_by_animal_behavior_y=value.png"))
plt.savefig(os.path.join(plotfolder, "coefficients_boxplots_by_animal_behavior_y=value.svg"))

# ------------------------------------

# Second plot with y = "ValueAbs"
g = sns.catplot(x="Coefficient", y="ValueAbs", hue="Animal", col="Behavior", data=Cdf, kind="box", height=4, aspect=2, legend_out=False, col_wrap=3)
g.set_xticklabels(rotation=45)

for ax in g.axes:
    ax.set_ylim(0, 1)
    ax.axhline(0, color="black", linestyle="--")
    ax.axvline(11.5, color="black", linestyle="--")
    for bar in ax.patches:
        bar.set_edgecolor("none")

# Save plots
g.savefig(os.path.join(plotfolder, "coefficients_boxplots_by_animal_behavior_y=abs.png"))
g.savefig(os.path.join(plotfolder, "coefficients_boxplots_by_animal_behavior_y=abs.svg"))

# ------------------------------------

# Line plot with y = "ValueAbs"
g = sns.relplot(x="Coefficient", y="ValueAbs", col="Behavior", data=Cdf, kind="line", height=4, aspect=2, col_wrap=3, hue="orthostate", ci="sd", facet_kws=dict(sharey=False))
g.set_xticklabels(rotation=45)
# Additional configurations
for ax in g.axes:
    ax.axhline(0, color="black", linestyle="--")
    for bar in ax.patches:
        bar.set_edgecolor("none")
for ax in g.axes:
    plt.sca(ax)
    mplcyberpunk.add_glow_effects()
    # For each subplot, let's find the maximum y-value for the lines and
    # draw an horizontal line at that value
    # ymax = max([line.get_ydata().max() for line in ax.lines])
    # ax.axhline(ymax, color="black", linestyle="--")
g.add_legend()
g.fig.savefig(os.path.join(plotfolder, "coefficients_lineplots_by_behavior_y=abs.png"))
g.fig.savefig(os.path.join(plotfolder, "coefficients_lineplots_by_behavior_y=abs.svg"))
g.fig.savefig(os.path.join(plotfolder, "coefficients_lineplots_by_behavior_y=abs.pdf"))

# ------------------------------------

# Line plot with y = "ValueAbs"
g = sns.relplot(x="Coefficient", y="ValueAbs", row="Animal", col="Behavior", data=Cdf, kind="line", height=4, aspect=2, hue="orthostate", ci="sd", facet_kws=dict(sharey=False))
g.set_xticklabels(rotation=45)
# Additional configurations
for ax in g.axes:
    ax.axhline(0, color="black", linestyle="--")
    for bar in ax.patches:
        bar.set_edgecolor("none")
for ax in g.axes.flat:
    plt.sca(ax)
    mplcyberpunk.add_glow_effects()
    # For each subplot, let's find the maximum y-value for the lines and
    # draw an horizontal line at that value
    # ymax = max([line.get_ydata().max() for line in ax.lines])
    # ax.axhline(ymax, color="black", linestyle="--")
g.add_legend()
# Format the title strings to have "\n" instead of "|" between row and column
for ax in g.axes.flat:
    ax.set_title(ax.get_title().replace("|", "\n"))
    # Switch order of row and column in title
    ax.set_title(ax.get_title().split("\n")[1] + "\n" + ax.get_title().split("\n")[0])
plt.subplots_adjust(top=0.95)
g.fig.savefig(os.path.join(plotfolder, "coefficients_lineplots_by_animal_behavior_y=abs.png"))
g.fig.savefig(os.path.join(plotfolder, "coefficients_lineplots_by_animal_behavior_y=abs.svg"))
g.fig.savefig(os.path.join(plotfolder, "coefficients_lineplots_by_animal_behavior_y=abs.pdf"))

# Let's get the ratio of the aligned versus orthogonal coefficients
# for each animal and behavior

Cdfp = Cdf.pivot_table(index=["Animal", "Behavior"], columns="orthostate", values="ValueAbs", aggfunc="mean")
Cdfp.loc[:,"ratio"] = Cdfp.Aligned / Cdfp.Orthogonal
Cdfp.reset_index(inplace=True)

sns.relplot(x="Behavior", y="ratio", data=Cdfp, kind="line", height=4, aspect=2, ci="sd")
plt.axhline(1.0, color="black", linestyle="--")
mplcyberpunk.add_glow_effects(gradient_fill=True)
plt.savefig(os.path.join(plotfolder, "ratio_aligned_orthogonal.png"))
plt.savefig(os.path.join(plotfolder, "ratio_aligned_orthogonal.svg"))
plt.savefig(os.path.join(plotfolder, "ratio_aligned_orthogonal.pdf"))


sns.catplot(x="ratio", data=Cdfp, kind="box", height=4, aspect=2, ci="sd")
mplcyberpunk.add_glow_effects(gradient_fill=True)
plt.savefig(os.path.join(plotfolder, "ratio_aligned_orthogonal_all_behaviors.png"))
plt.savefig(os.path.join(plotfolder, "ratio_aligned_orthogonal_all_behaviors.svg"))
plt.savefig(os.path.join(plotfolder, "ratio_aligned_orthogonal_all_behaviors.pdf"))
plt.suptitle("Ratio of Aligned versus Orthogonal Coefficients")
plt.axvline(x=1.0, color="black", linestyle="--", linewidth=3)

g=sns.catplot(x="ratio", data=Cdfp, kind="box", height=4, aspect=2, ci="sd", row="Animal")
plt.suptitle("Ratio of Aligned versus Orthogonal Coefficients")
g.map(lambda **kws: plt.axvline(x=1.0, color="black", linestyle="--", linewidth=3))
plt.subplots_adjust(top=0.90)
plt.savefig(os.path.join(plotfolder, "ratio_aligned_orthogonal_all_behaviors_by_animal.png"))
plt.savefig(os.path.join(plotfolder, "ratio_aligned_orthogonal_all_behaviors_by_animal.svg"))
plt.savefig(os.path.join(plotfolder, "ratio_aligned_orthogonal_all_behaviors_by_animal.pdf"))
plt.close()

g=sns.catplot(x="ratio", data=Cdfp, kind="box", height=4, aspect=2, ci="sd", row="Behavior")
plt.suptitle("Ratio of Aligned versus Orthogonal Coefficients")
g.map(lambda **kws: plt.axvline(x=1.0, color="black", linestyle="--", linewidth=3))
plt.subplots_adjust(top=0.90)
plt.savefig(os.path.join(plotfolder, "ratio_aligned_orthogonal_all_animals_by_behavior.png"))
plt.savefig(os.path.join(plotfolder, "ratio_aligned_orthogonal_all_animals_by_behavior.svg"))
plt.savefig(os.path.join(plotfolder, "ratio_aligned_orthogonal_all_animals_by_behavior.pdf"))

