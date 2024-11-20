import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pylab as plt
from matplotlib.backends.backend_pdf import PdfPages

# from Arin
# ---
# protein_fc:  fold change in enzyme usage fluxes during pure protein
# synthesis from the high-glucose condition to the parallel-favouring
# condition
#
# how enzyme usage fluxes change from pure carbohydrate synthesis to
# pure protein synthesis:
# C rich: highglucose_carb2prot_fc
# C poor: parallel_carb2prot_fc

df = pd.read_csv("fluxes.tsv", sep="\t")


def define_entry(e):
    """Find out what an entry does."""
    res = df[df["Entry Name"] == f"{e}_YEAST"]["Protein names"].values[0]
    for res_split in res.split(","):
        print(res_split.strip())


# tidy
sdf = df.copy()
sdf["entry_name"] = sdf["Entry Name"]
sdf.drop(["UNIPROT ID", "Protein names", "Entry Name"], axis=1, inplace=True)
sdf["entry_name"] = sdf["entry_name"].str.split("_").str[0]
sdf = sdf.rename(columns={"Gene Names": "gene_names"})
col = sdf.pop("entry_name")
sdf.insert(0, "entry_name", col)

# convert to long format
mydict = {
    "entry": np.tile(sdf.entry_name.values, 4),
    "flux": np.hstack(
        (
            sdf.carbohydrate_fluxes_highglucose.values,
            sdf.protein_fluxes_highglucose.values,
            sdf.carbohydrate_fluxes_parallel.values,
            sdf.protein_fluxes_parallel.values,
        )
    ),
    "type": (["carb"] * sdf.shape[0] + ["prot"] * sdf.shape[0]) * 2,
    "condition": ["rich"] * sdf.shape[0] * 2 + ["poor"] * sdf.shape[0] * 2,
}
ndf = pd.DataFrame.from_dict(mydict)
# sort alphabetically and numerically
ndf = ndf.sort_values(
    "entry",
    key=lambda x: x.str.extract("(\D+)", expand=False).str.lower()
    + x.str.extract("(\d+)", expand=False).fillna("0").str.zfill(10),
)

# define list of electron transport chain enzymes
with open("etc_enzymes.txt") as fobj:
    etc_enzymes = [line.rstrip() for line in fobj]

# bar plots
for c in ["rich", "poor"]:
    # choose fluxes relevant to nutrient condition
    ndf_condition = ndf[ndf.condition == c]
    # get list of enzymes
    enzymes = ndf_condition["entry"].unique().tolist()
    # move electron transport enzymes to top
    etc = [x for x in enzymes if x in etc_enzymes]
    non_etc = [x for x in enzymes if x not in etc_enzymes]
    enzymes = etc + non_etc

    # draw
    plt.figure(figsize=(5, 10))
    # bars
    sns.barplot(ndf_condition, y="entry", x="flux", hue="type", order=enzymes)
    # highlight background of electron transport chain enzymes
    plt.axhspan(-0.5, len(etc_enzymes) - 0.5, facecolor="#c3d2a8", zorder=-100)
    plt.title(c)
    plt.show(block=False)

# Print drawings
pdf_filename = f"fluxes.pdf"
with PdfPages(pdf_filename) as pdf:
    for fig in range(1, plt.gcf().number + 1):
        pdf.savefig(fig)
# Close all figures
plt.close("all")
