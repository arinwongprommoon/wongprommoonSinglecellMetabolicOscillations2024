import colorcet as cc
import matplotlib.pylab as plt
from wela.dataloader import dataloader
from wela.plotting import plot_cuml_divisions_per_cell

from data_dict import data_dict

# FY4
i = 0
plt.figure(figsize=(10, 5))
for data_id in data_dict:
    if "19972" in data_dict[data_id][0] or "20212" in data_dict[data_id][0]:
        continue
    elif "613" in data_dict[data_id][0]:
        data_dict[data_id][2] = "2% glu K"
    dl = dataloader(wdir="tsv_data/individual", ls=False)
    dl.load(data_dict[data_id][0], use_tsv=True)
    for group in dl.df.group.unique():
        if group == "htb2mCherry":
            t, buddings = dl.get_time_series("buddings", group=group)
            plot_cuml_divisions_per_cell(
                t,
                buddings,
                label=f"{data_dict[data_id][2]}: {group}",
                col=cc.glasbey_bw[i],
            )
            i += 1
plt.xlabel("time")
plt.ylabel("mean number of cumulative divisions")
plt.legend(bbox_to_anchor=(1.38, 0.9))
plt.tight_layout()
plt.show(block=False)

# BY and mutants
i = 0
plt.figure(figsize=(10, 5))
for data_id in data_dict:
    if "19972" in data_dict[data_id][0] or "20212" in data_dict[data_id][0]:
        continue
    elif "613" in data_dict[data_id][0]:
        data_dict[data_id][2] = "2% glu K"
    dl = dataloader(wdir="tsv_data/individual", ls=False)
    dl.load(data_dict[data_id][0], use_tsv=True)
    for group in dl.df.group.unique():
        if group != "htb2mCherry":
            t, buddings = dl.get_time_series("buddings", group=group)
            plot_cuml_divisions_per_cell(
                t,
                buddings,
                label=f"{data_dict[data_id][2]}: {group}",
                col=cc.glasbey_bw[i],
            )
            i += 1
plt.xlabel("time")
plt.ylabel("mean number of cumulative divisions")
plt.legend(bbox_to_anchor=(1.35, 0.9))
plt.tight_layout()
plt.show(block=False)
