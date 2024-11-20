This repository contains the code used for the manuscript, "Single-cell metabolic oscillations are pervasive and may alleviate a proteome constraint", authored by Arin Wongprommoon, Alán F. Muñoz, Diego A. Oyarzún, and Peter S. Swain.

# Directory structure

This repository contains three directories:

- "timeseries" : Used for analysing and visualising time series.
- "phase" : Used for analysing phases of time series.
- "fba" : Used for flux balance analysis.

# General requirements

The code has been tested on Python 3.8.

# timeseries

## Requirements

This directory requires these Python dependencies:

- `aliby` and its dependencies : documented on [https://pypi.org/project/aliby/], source code available at [https://gitlab.com/aliby/aliby]
- `matplotlib`
- `numpy`
- `wela` and its dependencies : documented on [https://git.ecdf.ed.ac.uk/swain-lab/wela]

Users will also need time-series data and classifications of time-series data from [https://doi.org/10.7488/ds/7840].

## Usage

The `strain_report.py` script from this directory can be adapted to:

- draw a plot showing the fluorescence (flavin or mCherry) of a single cell over time, annotated with birth events (relevant to Figures 1A and S2C), or
- draw a heatmap showing how the flavin fluorescence of multiple cells in an experiment changes over time (relevant to Figure 1E).

The rest of the files are modules containing functions and variables used by the main plotting script.

# phase

## Requirements

This directory requires these Python dependencies:

- `astropy`
- `colorcet`
- `matplotlib`
- `numpy`
- `pandas`
- `scipy`
- `seaborn`
- `wela` and its dependencies : documented on [https://git.ecdf.ed.ac.uk/swain-lab/wela]

Users will also need time-series data and classifications of time-series data from [https://doi.org/10.7488/ds/7840].

## Usage

This directory produces plots relevant to Figures 1B-D, 2, S1, S2A-B, and S3-6.

To draw all plots, run `python plot_phase.py`.  This script may be time-consuming -- allow approximately 20 minutes.  This script will produce a PDF file called `all_plots.pdf`.

`data_dict.py` defines the list of datasets (time-series data and classification of time-series data) to be used for time series analysis and drawing plots.  For the code in this directory to work, these datasets must be in a "`tsv_data`" directory in the same working directory as the Python scripts.

The rest of the files are modules containing functions and variables used by the main plotting script.

# fba

## Requirements

This directory requires these Python dependencies:

- `cobra`
- `matplotlib`
- `numpy`
- `seaborn`
- `wrapt_timeout_decorator`

Users will also need the ecYeast8.6.0 model, available at  [https://github.com/SysBioChalmers/ecModels/blob/update/ecYeastGEM/v8.6.0/2022-07-03-20-48-59/ecYeastGEM/model/ecYeastGEM_batch.xml].  Save this file as `ecYeastGEM_batch_8-6-0.xml` in the same directory as the Python scripts for the scripts to work.

Users will also need flux data from [https://doi.org/10.7488/ds/7840], saved as `fluxes.tsv`.

## Usage

This directory produces plots relevant to Figures 3B-C and Figure 4.

To draw the bar plots in Figure 3B, run `python plot_strains.py`.

To draw the heatmap in Figure 3C, run `python plot_enzyme_usage.py`.

To draw the heatmaps in Figure 4A, first run `python grid_ablation.py` to run metabolic simulations and produce output *.pkl files.  This step may take several hours, or overnight.  Then run `python plot_grid_ablation.py` to draw the heatmaps.  Modify the `model_options` variable to switch between glucose and pyruvate as the carbon source.

To draw the bar plots in Figure 4B, run `python plot_fluxes.py`.  This requires the `fluxes.tsv` from the DataStore and the list of electron transport chain enzymes, `etc_enzymes.txt`.

To obtain respiration and fermentation fluxes relevant to the section, "A constrained proteome usually favours sequential biosynthesis", run `python respiration_fermentation.py`.  The output file is `respiration_and_fermentation.csv`, with the columns:

- "c_source" : Carbon source, "glc" for glucose, "pyr" for pyruvate.
- "condition" : Nutrient condition defined by carbon and nitrogen uptake, "seqbetter" for the condition favouring sequential synthesis, "parbetter" for the condition favouring parallel synthesis.
- "biomass_component" : The simulation simulates the cell producing this particular biomass component. 
- "exchange_species" : Chemical species that defines the exchange rate, e.g. "O2" denotes oxygen exchange rate (uptake).
- "flux" : Flux of the exchange reaction defined by "exchange_species".

The rest of the files are modules containing functions and variables used by the plotting scripts.
