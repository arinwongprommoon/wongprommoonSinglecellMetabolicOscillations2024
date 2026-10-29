This repository contains the code used for the manuscript, "Single-cell metabolic oscillations may alleviate a proteome constraint", authored by Arin Wongprommoon, Diego A. Oyarzún, and Peter S. Swain.

# Directory structure

The code to reproduce figures is in the "fba" directory.

# Requirements

The code has been tested on Python 3.8.

This directory requires these Python dependencies:

- `cobra`
- `matplotlib`
- `numpy`
- `seaborn`
- `wrapt_timeout_decorator`

Users will also need the ecYeast8.6.0 model, available at  [https://github.com/SysBioChalmers/ecModels/blob/update/ecYeastGEM/v8.6.0/2022-07-03-20-48-59/ecYeastGEM/model/ecYeastGEM_batch.xml].  Save this file as `ecYeastGEM_batch_8-6-0.xml` in the same directory as the Python scripts for the scripts to work.

Users will also need flux data from [https://doi.org/10.7488/ds/7845], saved as `fluxes.tsv`.  This file can also be generated using `script_compare_enz_use.py`, which operates on the UNIPROT-derived lookup table `idmapping_2024_07_27.tsv`.

## Usage

This directory produces plots relevant to Figures 1B-C and Figures 2A-2B.

To draw the bar plots in Figure 1B, run `python plot_strains.py`.

To draw the heatmap in Figure 1C, run `python plot_enzyme_usage.py`.

To draw the heatmaps in Figure 2A, first run `python grid_ablation.py` to run metabolic simulations and produce output *.pkl files.  This step may take several hours, or overnight, so the pkl files are included in this repo.  Then run `python plot_grid_ablation.py` to draw the heatmaps.  Modify the `model_options` variable to switch between glucose and pyruvate as the carbon source.

To draw the bar plots in Figure 2B, run `python plot_fluxes.py`.  This requires the `fluxes.tsv` from the DataStore and the list of electron transport chain enzymes, `etc_enzymes.txt`.

To obtain respiration and fermentation fluxes relevant to the section, "A constrained proteome usually favours sequential biosynthesis", run `python respiration_fermentation.py`.  The output file is `respiration_and_fermentation.csv`, with the columns:

- "c_source" : Carbon source, "glc" for glucose, "pyr" for pyruvate.
- "condition" : Nutrient condition defined by carbon and nitrogen uptake, "seqbetter" for the condition favouring sequential synthesis, "parbetter" for the condition favouring parallel synthesis.
- "biomass_component" : The simulation simulates the cell producing this particular biomass component. 
- "exchange_species" : Chemical species that defines the exchange rate, e.g. "O2" denotes oxygen exchange rate (uptake).
- "flux" : Flux of the exchange reaction defined by "exchange_species".

The rest of the files are modules containing functions and variables used by the plotting scripts.
