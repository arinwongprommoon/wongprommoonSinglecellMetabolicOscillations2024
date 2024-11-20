import matplotlib.pylab as plt
from wela.dataloader import dataloader
from wela.fft_filtering import flavin_filter
from wela.plotting import bud_to_bud_plot

data_dict = {
    0: [
        "26643_2022_05_23_flavin_htb2_glucose_20gpL_01_00_htb2mCherry.tsv",
        "2% glu",
        ["htb2mCherry"],
    ],
    1: [
        "31594_2022_09_25_flavin_htb2_pyruvate_20gpL_01_03_htb2mCherry.tsv",
        "2% pyr",
        ["htb2mCherry"],
    ],
    2: [
        "1649_2023_04_18_flavin_by4742swain_by4742morgan_tsa1tsa2morgan_lysmedia_02_00.tsv",
        "1% glu",
        ["by4742morgan", "tsa1tsa2morgan"],
    ],
    3: [
        "1253_2023_03_23_flavin_by4742swain_by4742morgan_tsa1tsa2morgan_lysmedia_01_00.tsv",
        "1% glu",
        ["by4742morgan", "tsa1tsa2morgan"],
    ],
    4: [
        "20016_2021_07_09_flavin_by4741_zwf1egf_01_00.tsv",
        "1% glu",
        ["by4741", "zwf1egf"],
    ],
    5: [
        "613_2023_01_30_flavin_fy4_htb2_k_deficient_01_01_htb2mCherry.tsv",
        "2% glu K -> 2% glu no K @ 6 h -> 2% glu K @ 16 h",
        ["htb2mCherry"],
    ],
}

times = [4, 9, 12]


# plot time series at specified time points
for data_id in data_dict:
    dl = dataloader(wdir="tsv_data/individual", ls=False)
    dl.load(data_dict[data_id][0], use_tsv=True)
    odf = dl.df
    for t in times:
        title = dl.dataname.split("_")[0] + " - " + data_dict[data_id][1] + f" - t={t} "
        if data_dict[data_id][2][0] == "htb2mCherry":
            # htb2
            fig, ax1 = plt.subplots()
            plt.sca(ax1)
            bud_to_bud_plot(
                t,
                "flavin",
                dl,
                group=data_dict[data_id][2][0],
                colour="b",
                df=odf,
                filter_func=flavin_filter,
                show_figure=False,
            )
            ax2 = ax1.twinx()
            plt.sca(ax2)
            bud_to_bud_plot(
                t,
                "mCherry",
                dl,
                group=data_dict[data_id][2][0],
                colour="r",
                df=odf,
                filter_func=None,
                show_figure=False,
            )
        else:
            # mutants
            plt.figure()
            # wild type
            bud_to_bud_plot(
                t,
                "flavin",
                dl,
                group=data_dict[data_id][2][0],
                colour="b",
                df=odf,
                filter_func=flavin_filter,
                show_figure=False,
            )
            # mutant
            bud_to_bud_plot(
                t,
                "flavin",
                dl,
                group=data_dict[data_id][2][1],
                colour="m",
                df=odf,
                filter_func=flavin_filter,
                show_figure=False,
            )
            title += f" {data_dict[data_id][2][1]}"
        plt.title(title)
        plt.show(block=False)
