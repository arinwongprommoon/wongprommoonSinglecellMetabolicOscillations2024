import matplotlib
from wela.dataloader import dataloader
from wela.fft_filtering import flavin_filter
from wela.plotting import kymograph
from wela.sorting import sort_by_budding

from data_dict import data_dict, get_oscillators

use_oscillators = True
data_id = 1

dl = dataloader(wdir="tsv_data/individual", ls=False)
dl.load(data_dict[data_id][0], use_tsv=True)
osc = get_oscillators(data_dict[data_id][1])
title = data_dict[data_id][2]
group = "htb2mCherry"

if use_oscillators:
    odf = dl.df[dl.df.id.isin(osc[group][0])]
else:
    odf = dl.df
_, buddings = dl.get_time_series("buddings", df=odf)
t, data = dl.get_time_series("flavin", df=odf)


sort_order = sort_by_budding(buddings)

kymograph(
    odf,
    group=group,
    hue="flavin",
    title=title,
    cmap=matplotlib.cm.RdBu,
    filterfunc=flavin_filter,
    standardscale=True,
    robust=True,
    buddings=True,
    figsize=(12, 6),
    sort_order=sort_order,
)
