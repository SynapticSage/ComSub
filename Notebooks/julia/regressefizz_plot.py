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
plotfolder = os.path.join(plotfolder, "julia_regressefizz")
if not os.path.exists(plotfolder):
    os.makedirs(plotfolder)

# Read the CSV file
Cdf = pd.read_csv(os.path.join(folder, "efizz_coefs_df.csv"))
assert((Cdf.Coefficient2 != "avg").all(), "Coefficient2 should not be avg")

# Filter out rows where Coefficient is "Intercept"
# Cdf.query('target_prop1 != "Intercept" & (orthostate == "Orthogonal" | orthostate == "Aligned")', inplace=True)

# ------------------------------------


# Box plot of coefficients by animal and behavior
g1 = sns.relplot(x="Coefficient", y="ValueAbs", hue="Coefficient1",
                 row="target_prop1", data=Cdf, kind="line", height=4, aspect=2,
                 col="target_prop2")
g1.set_xticklabels(rotation=45)
# Change title to \n instead of | between title pieces
for ax in g1.axes.flat:
    ax.set_title(ax.get_title().replace('|', '\n'))
# Saving the plot
# plt.savefig("coefficients_boxplots_by_animal_behavior.png")

from tqdm import tqdm
for coef2 in tqdm(Cdf.Coefficient2.unique()):
    cdf = Cdf.query('Coefficient2 == @coef2')
    assert len(cdf) > 0
    g1 = sns.catplot(x="Coefficient", y="ValueAbs", hue="Coefficient1",
                     row="target_prop1", data=cdf, kind="box", height=4, aspect=2,
                     col="target_prop2", dodge=False)
    g1.set_xticklabels(rotation=45)
    # Change title to \n instead of | between title pieces
    for ax in g1.axes.flat:
        ax.set_title(ax.get_title().replace('|', '\n'))
        # decrease thickness of box lines
        for box in ax.artists:
            box.set_linewidth(0.5)
    g1.fig.suptitle("Coefficient2 = " + coef2)
    g1.fig.subplots_adjust(right=0.8, top=0.9)
    g1.fig.show()
    g1.add_legend()
    # Saving the plot
    plt.savefig(os.path.join(plotfolder, "coefficients_boxplots_by_animal_behavior_Coefficient2=" + coef2 + ".png"))
    plt.savefig(os.path.join(plotfolder, "coefficients_boxplots_by_animal_behavior_Coefficient2=" + coef2 + ".pdf"))

from tqdm import tqdm
for coef1 in tqdm(Cdf.Coefficient1.unique()):
    cdf = Cdf.query('Coefficient1 == @coef1')
    assert len(cdf) > 0
    g1 = sns.catplot(x="Coefficient", y="ValueAbs", hue="Coefficient2",
                     row="target_prop1", data=cdf, kind="box", height=4, aspect=2,
                     col="target_prop2", dodge=False)
    g1.set_xticklabels(rotation=45)
    # Change title to \n instead of | between title pieces
    for ax in g1.axes.flat:
        ax.set_title(ax.get_title().replace('|', '\n'))
        # decrease thickness of box lines
        for box in ax.artists:
            box.set_linewidth(0.5)
    g1.fig.suptitle("Coefficient1 = " + coef1)
    g1.fig.subplots_adjust(right=0.8, top=0.9)
    g1.fig.show()
    g1.add_legend()
    # Saving the plot
    plt.savefig(os.path.join(plotfolder, "coefficients_boxplots_by_animal_behavior_Coefficient1=" + coef1 + ".png"))
    plt.savefig(os.path.join(plotfolder, "coefficients_boxplots_by_animal_behavior_Coefficient1=" + coef1 + ".pdf"))


# Pivot Coefficient1 
cdfp = Cdf.pivot_table(columns="target_prop1", index=["Coefficient", "Coefficient2", "Coefficient1", "target_prop2"], values="ValueAbs")
cdfp.reset_index(inplace=True)
cdfp.loc[:,"ratio"] = cdfp.loc[:,"R"] / cdfp.loc[:,"Ru"]

# Box plot of coefficients by animal and behavior
g1 = sns.relplot(x="Coefficient", y="ratio", hue="Coefficient1",
                 row="target_prop2", data=cdfp, kind="line", height=4,
                 aspect=2)
