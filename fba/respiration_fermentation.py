#!/usr/bin/env python
import pandas as pd
from yeast8model import Yeast8Model

# # Construct models of cells of interest, optimise

glc_exch_rate = 16.89

wt = Yeast8Model("ecYeastGEM_batch_8-6-0.xml")
wt.model.reactions.get_by_id("r_1714").bounds = (-glc_exch_rate, 0)
wt.model.reactions.get_by_id("r_1714_REV").bounds = (0, glc_exch_rate)
wt.solution = wt.optimize()

# # Ablate
# ## Wild type

wt.ablation_result = wt.ablate()
wt_ablation_result_orig = wt.ablation_result.copy()


# ## Respiration and fermentation

# Vary nutrient sources

preset_models = {
    "glc_seqbetter": {
        "glc_exch_rate": 16.89,
        "pyr_exch_rate": None,
        "amm_exch_rate": None,
    },
    "glc_parbetter": {
        "glc_exch_rate": 0.194 * 8.6869,
        "pyr_exch_rate": None,
        "amm_exch_rate": 0.71 * 1.4848,
    },
    "pyr_seqbetter": {
        "glc_exch_rate": 0,
        "pyr_exch_rate": 8.8888,
        "amm_exch_rate": None,
    },
    "pyr_parbetter": {
        "glc_exch_rate": 0,
        "pyr_exch_rate": 0.839 * 4.4444,
        "amm_exch_rate": 0.903 * 1.0,
    },
}

biomass_components = ["carbohydrate", "protein"]

exchanges = {
    "O2": "r_1992_REV",
    "CO2": "r_1672",
    "ethanol": "r_1761",
}


def get_exchanges(model_name, model_options):
    # Set nutrient uptake rates based on conditions
    if model_options["glc_exch_rate"] is None:
        wt.model.reactions.get_by_id("r_1714").bounds = (-16.89, 0)
        wt.model.reactions.get_by_id("r_1714_REV").bounds = (0, 16.89)
    else:
        wt.model.reactions.get_by_id("r_1714").bounds = (
            -model_options["glc_exch_rate"],
            0,
        )
        wt.model.reactions.get_by_id("r_1714_REV").bounds = (
            0,
            model_options["glc_exch_rate"],
        )

    if model_options["pyr_exch_rate"] is None:
        pass
    else:
        # wt.model.reactions.get_by_id("r_1714").bounds = (0, 0)
        wt.model.reactions.get_by_id("r_2033").bounds = (
            -model_options["pyr_exch_rate"],
            0,
        )
        wt.model.reactions.get_by_id("r_2033_REV").bounds = (
            0,
            model_options["pyr_exch_rate"],
        )

    if model_options["amm_exch_rate"] is None:
        pass
    else:
        wt.model.reactions.get_by_id("r_1654").bounds = (
            -model_options["amm_exch_rate"],
            0,
        )
        wt.model.reactions.get_by_id("r_1654_REV").bounds = (
            0,
            model_options["amm_exch_rate"],
        )

    # Simulate
    wt.solution = wt.optimize()
    wt.ablation_result = wt.ablate()

    # Store data on fluxes
    c_source_list = []
    condition_list = []
    biomass_component_list = []
    exchange_species_list = []
    flux_list = []

    c_source, condition = model_name.split("_")

    for biomass_component in biomass_components:
        for exchange_species, exchange_id in exchanges.items():
            flux = wt.ablation_reaction_fluxes[biomass_component][exchange_id]

            c_source_list.append(c_source)
            condition_list.append(condition)
            biomass_component_list.append(biomass_component)
            exchange_species_list.append(exchange_species)
            flux_list.append(flux)

    df = pd.DataFrame(
        {
            "c_source": c_source_list,
            "condition": condition_list,
            "biomass_component": biomass_component_list,
            "exchange_species": exchange_species_list,
            "flux": flux_list,
        }
    )
    return df


df_list = []
for key, val in preset_models.items():
    wt = Yeast8Model("ecYeastGEM_batch_8-6-0.xml")
    df = get_exchanges(key, val)
    df_list.append(df)

res_df = pd.concat(df_list)
res_df.to_csv("respiration_and_fermentation.csv")
