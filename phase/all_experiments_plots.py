import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set_style("whitegrid")


def augment_data(mean, std):
    """Augment data so seaborn generates appropriate error bars."""
    if np.isnan(mean) or np.isnan(std):
        return [np.nan, np.nan, np.nan]
    else:
        data = [
            mean,
            mean + np.sqrt(3 / 2) * std,
            mean - np.sqrt(3 / 2) * std,
        ]
        return data


def create_full_df(all_phases_dfs):
    """Create a dataframe for the single-cell periods and ratios."""
    full_df = pd.concat([df for df in all_phases_dfs], ignore_index=True)
    # drop phase and all duplicates of periods
    full_df = full_df.drop(
        columns=[col for col in full_df.columns if "phase" in col]
    ).drop_duplicates()
    # drop IDs because their use was to find duplicates
    full_df = full_df.drop(columns=[col for col in full_df.columns if "_id" in col])
    # add ratios
    full_df["flavin_to_buddings_period_ratio"] = np.log2(
        full_df.flavin_period / full_df.buddings_period
    )
    full_df["htb2_to_buddings_period_ratio"] = np.log2(
        full_df.htb2_period / full_df.buddings_period
    )
    full_df["htb2_to_flavin_period_ratio"] = np.log2(
        full_df.htb2_period / full_df.flavin_period
    )
    # rearrange for plotting
    sdf = pd.melt(
        full_df,
        id_vars=["condition", "group", "omero_no"],
        var_name="type",
        value_name="value",
    )
    sdf[["variable", "stat"]] = sdf.type.str.rsplit("_", n=1, expand=True)
    sdf = sdf.drop(columns=["type"])
    return sdf


def plot_periods_ratios(summary_df):
    """Plot periods and ratios of periods from summary_df."""
    for xv, hv in zip(["period", "ratio"], ["period_variable", "ratio_type"]):
        # wild type
        plt.figure()
        ax = sns.barplot(
            summary_df[summary_df.group == "htb2mCherry"],
            x=xv,
            y="condition",
            hue=hv,
            errorbar="sd",
        )
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
        ax.grid()
        plt.tight_layout()
        plt.show(block=False)
        # mutants
        for yv in ["group", "id"]:
            plt.figure()
            ax = sns.barplot(
                summary_df[
                    (summary_df.condition == "1% glu")
                    & summary_df.period_variable.str.contains("flavin|buddings")
                ],
                x=xv,
                y=yv,
                hue=hv,
                errorbar="sd",
            )
            sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
            ax.grid()
            plt.tight_layout()
            plt.show(block=False)


def plot_single_cell_data(all_phases_dfs):
    """Plot single-cell periods and ratios of periods."""
    # combine single-cell data from all experiments
    full_df = create_full_df(all_phases_dfs)
    for stat in ["period", "ratio"]:
        # for wild type
        sdf = full_df[(full_df.group == "htb2mCherry") & (full_df.stat == stat)]
        g = sns.catplot(
            data=sdf,
            x="value",
            y="condition",
            hue="variable",
            kind="box",
            showfliers=False,
        )
        g.fig.suptitle(stat)
        plt.show(block=False)
        # distributions
        g = sns.displot(
            data=sdf,
            x="value",
            hue="variable",
            col="condition",
            height=4,
            aspect=0.7,
            kind="kde",
        )
        g.fig.suptitle(stat)
        plt.show(block=False)
        # for mutants
        sdf = full_df[
            (full_df.condition == "1% glu")
            & full_df.variable.str.contains("flavin|buddings")
            & (full_df.stat == stat)
        ]
        g = sns.catplot(
            data=sdf,
            x="value",
            y="group",
            hue="variable",
            kind="box",
            showfliers=False,
        )
        g.fig.suptitle(stat)
        plt.show(block=False)
        # distributions
        g = sns.displot(
            data=sdf,
            x="value",
            hue="variable",
            col="group",
            height=4,
            aspect=0.7,
            kind="kde",
        )
        g.fig.suptitle(stat)
        plt.show(block=False)


def catplot_from_period_dict(period_dict):
    """Plot using most likely periods."""
    df_dict = {
        key: []
        for key in [
            "omero_no",
            "condition",
            "group",
            "period",
            "replicate",
            "type",
        ]
    }
    for key in period_dict:
        if key.endswith("_period"):
            for i in range(len(period_dict[key])):
                adata = augment_data(period_dict[key][i], period_dict[f"{key}_std"][i])
                for j, adata_pt in enumerate(adata):
                    df_dict["period"].append(adata_pt)
                    df_dict["replicate"].append(j)
                    for field in ["omero_no", "condition", "group"]:
                        df_dict[field].append(period_dict[field][i])
                    df_dict["type"].append(key.split("_")[0])
    df = pd.DataFrame.from_dict(df_dict)
    # data with htb2
    g = sns.catplot(
        data=df[df.group == "htb2mCherry"],
        kind="bar",
        y="condition",
        x="period",
        hue="type",
        errorbar="se",
    )
    g.fig.suptitle("htb2mCherry")
    plt.show(block=False)
    # data without htb2
    plt.show(block=False)
    g = sns.catplot(
        data=df[df.condition == "1% glu"],
        kind="bar",
        y="group",
        x="period",
        hue="type",
        errorbar="se",
    )
    g.fig.suptitle("1% glu")
    plt.show(block=False)


def plot_from_period_dict(period_dict):
    """Plot using most likely periods."""
    df = pd.DataFrame.from_dict(period_dict)
    # skip pyruvate to spread out data
    df = df[df.condition != "2% pyr"]
    xv = "flavin_period"
    yv = "buddings_period"
    fig, ax = plt.subplots()
    sns.scatterplot(data=df, x=xv, y=yv, ax=ax, hue="group", style="condition", s=100)
    ax.errorbar(
        df[xv],
        df[yv],
        xerr=df[f"{xv}_std"],
        yerr=df[f"{yv}_std"],
        fmt="none",
        ecolor="gray",
        alpha=0.2,
        capsize=5,
    )
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    min_val = min(xlim[0], ylim[0])
    max_val = max(xlim[1], ylim[1])
    ax.set_xlim(min_val, max_val)
    ax.set_ylim(min_val, max_val)
    ax.set_aspect("equal")
    plt.tight_layout()
    plt.show(block=False)
