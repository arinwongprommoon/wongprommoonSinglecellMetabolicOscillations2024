#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

prop_cycle = plt.rcParams["axes.prop_cycle"]
default_mpl_colors = prop_cycle.by_key()["color"]


def heatmap_ablation_grid(
    ax,
    exch_rate_dict,
    ratio_array,
    largest_component_array=None,
    percent_saturation=False,
    saturation_point=(None, None),
    saturation_grid=False,
    vmin=0,
    vmax=2,
    showticklabels=True,
    center=None,
    cmap="RdBu_r",
    cbar_label="ratio",
):
    """Draw heatmap from 2d ablation grid

    Parameters
    ----------
    ax : matplotlib.pyplot.Axes object
        Axes to draw heatmap on.
    exch_rate_dict : dict
        dict that stores the two exchange reactions to vary and the uptake
        rate values to use.  It should be in this format:

        d = {
            'r_exch_rxn_1' : <array-like>,
            'r_exch_rxn_2' : <array-like>,
            }

    ratio_array : numpy.ndarray (2-dimensional)
        Array of ablation ratios, output from ablation_grid()
    largest_component_array : numpy.ndarray (2-dimensional), optional
        Array of largest biomass components, output from ablation_grid()
    percent_saturation : bool, optional
        Whether to scale axis labels so that the numbers displayed are percent
        of saturation.  Default False.
    saturation_point : (float, float) tuple, optional
        Values of exchange fluxes to use as saturation.  If either element is
        None, use the max as saturation.
    saturation_grid : bool, optional
        Whether to draw grid lines to show where saturation is.  Default False.
    vmin : float, optional
        Minimum of range for colour bar.  Default 0.
    vmax : float, optional
        Maximum of range for colour bar.  Default 2.
    showticklabels : bool, optional
        Whether to show x/y-axis tick labels.  Default True.
    cbar_label : string, optional
        Label for colour bar.  Default "ratio".

    Examples
    --------
    FIXME: Add docs.

    """
    # If largest_component_array is supplied, use it as text labels on heatmap.
    # This design takes advantage of seaborn.heatmap(annot=None) being default.
    if largest_component_array is None:
        annot_input = largest_component_array
    # TODO: Improve error-handling by checking if this is a 2D numpy array
    else:
        annot_input = np.rot90(largest_component_array)

    # Define x & y tick labels
    heatmap_xticklabels = list(exch_rate_dict.values())[0].copy()
    heatmap_yticklabels = list(exch_rate_dict.values())[1][::-1].copy()
    saturation_x = saturation_point[0]
    saturation_y = saturation_point[1]
    # ... depending on saturation-related arguments
    if percent_saturation:
        if saturation_x is not None:
            heatmap_xticklabels /= saturation_x
        else:
            heatmap_xticklabels /= np.max(heatmap_xticklabels)
        heatmap_xticklabels *= 100

        if saturation_y is not None:
            heatmap_yticklabels /= saturation_y
        else:
            heatmap_yticklabels /= np.max(heatmap_yticklabels)
        heatmap_yticklabels *= 100

        # and draw grid lines if specified
        # This only makes sense if percent_saturation is True.
        if saturation_grid:
            ax.axvline(
                np.searchsorted(heatmap_xticklabels, 100, side="left"),
                color="k",
                linewidth=1,
            )
            # doing this because y axis is defined 'in reverse' & to have line
            # position consistent with x axis
            ax.axhline(
                np.searchsorted(heatmap_yticklabels[::-1], 100, side="right"),
                color="k",
                linewidth=1,
            )

    # Draws heatmap.
    # Rounding directly on the x/yticklabels variables because of known
    # matplotlib-seaborn bug:
    # - https://github.com/mwaskom/seaborn/issues/1005
    # - https://stackoverflow.com/questions/63964006/round-decimal-places-seaborn-heatmap-labels
    # - https://stackoverflow.com/questions/50571592/matplotlib-formatstrformatter-returns-wrong-values
    sns.heatmap(
        data=np.rot90(ratio_array),
        annot=annot_input,
        xticklabels=np.around(heatmap_xticklabels, decimals=1),
        yticklabels=np.around(heatmap_yticklabels, decimals=1),
        robust=True,
        vmin=vmin,
        vmax=vmax,
        center=center,
        cmap=cmap,
        cbar_kws={"label": cbar_label},
        fmt="",
        ax=ax,
    )
    ax.set_xlabel(list(exch_rate_dict.keys())[0])
    ax.set_ylabel(list(exch_rate_dict.keys())[1])

    # Hide tick labels
    if not showticklabels:
        ax.set(xticklabels=[])
        ax.set(yticklabels=[])
