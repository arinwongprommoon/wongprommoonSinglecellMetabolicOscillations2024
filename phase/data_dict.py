import pandas as pd

data_dict = {
    0: [
        "26643_2022_05_23_flavin_htb2_glucose_20gpL_01_00_htb2mCherry.tsv",
        ["is26643_htb2mCherry_labels.csv"],
        "2% glu",
    ],
    1: [
        "31594_2022_09_25_flavin_htb2_pyruvate_20gpL_01_03_htb2mCherry.tsv",
        ["is31594_htb2mCherry_labels.csv"],
        "2% pyr",
    ],
    2: [
        "1649_2023_04_18_flavin_by4742swain_by4742morgan_tsa1tsa2morgan_lysmedia_02_00.tsv",
        ["st01649_tsa1tsa2morgan_labels.csv"],
        "1% glu",
    ],
    3: [
        "31492_2022_09_23_flavin_htb2_glucose_10mgpL_03_00_htb2mCherry.tsv",
        ["is31492_htb2mCherry_labels.csv"],
        "0.001% glu",
    ],
    4: [
        "1253_2023_03_23_flavin_by4742swain_by4742morgan_tsa1tsa2morgan_lysmedia_01_00.tsv",
        ["st01253_by4742swain_labels.csv", "st01253_tsa1tsa2morgan_labels.csv"],
        "1% glu",
    ],
    5: [
        "19972_2021_05_28_flavin_htb2_glucose_limitation_hard_Delft_01_01_htb2mCherry.tsv",
        ["is19972_htb2mCherry_labels.csv"],
        "0.75% glu -> 0% glu @ 7 h -> 0.75% glu @ 15 h",
    ],
    6: [
        "20016_2021_07_09_flavin_by4741_zwf1egf_01_00.tsv",
        ["is20016_by4741_labels.csv", "is20016_zwf1egf_labels.csv"],
        "1% glu",
    ],
    7: [
        "20212_2021_09_04_flavin_cenpkmata_tsa1tsa2_rim11_swe1_htb2_08_00.tsv",
        ["is20212_cenpkkoetter_labels.csv"],
        "CEN.PK 1% glu",
    ],
    8: [
        "27917_2022_06_14_flavin_htb2_glucose_20gpL_02_00_htb2mCherry.tsv",
        ["is27917_htb2mCherry_labels.csv"],
        "2% glu",
    ],
    9: [
        "613_2023_01_30_flavin_fy4_htb2_k_deficient_01_01_htb2mCherry.tsv",
        ["st00613_htb2mCherry_labels.csv"],
        "2% glu K -> 2% glu no K @ 6 h -> 2% glu K @ 16 h",
    ],
    10: [
        "27895_2022_06_10_flavin_htb2_glucose_10mgpL_02_00_htb2mCherry.tsv",
        ["is27895_htb2mCherry_labels.csv"],
        "0.001% glu",
    ],
    11: [
        "18617_2020_02_21_protAgg_downUpShift_2_0_2_pHluorin_Ura7HA_Ura7HR_00.tsv",
        ["is18617_pHluorin_labels.csv"],
        "2% glu -> 0% glu @ 5 h -> 2% glu @ 10 h",
    ],
}


def get_oscillators(datanames, direc="tsv_data/oscillating"):
    """Get list of oscillating cells."""
    res = {}
    for dataname in datanames:
        ldf = pd.read_csv(f"{direc}/{dataname}")
        ldf.insert(
            0,
            "id",
            ldf["position"]
            + ";"
            + ldf["cell_label"].astype(str)
            + ";"
            + ldf["trap"].astype(str),
        )
        oscillators = list(ldf.id[ldf.score == 1].values)
        nonoscillators = list(ldf.id[ldf.score == 0].values)
        res[dataname.split("_")[1]] = [oscillators, nonoscillators]
    return res
