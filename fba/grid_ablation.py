#!/usr/bin/env python
import numpy as np
import pickle
from yeast8model import Yeast8Model


# # Construct models of cells of interest
glc_exch_rate = 16.89
wt_ec = Yeast8Model("ecYeastGEM_batch_8-6-0.xml")
wt_ec.model.reactions.get_by_id("r_1714").bounds = (-glc_exch_rate, 0)
wt_ec.model.reactions.get_by_id("r_1714_REV").bounds = (0, glc_exch_rate)

# # Generate data

# `linspace` values are based on saturation exchange rates from the saturation
# curves.  These values may be different from the optimal uptake values from
# optimising the unmodified wild-type model.

# ### Glucose × ammonium
exch_rate_dict = {
    "r_1714": np.linspace(0, 2 * 8.6869, 32),  # glucose
    "r_1654": np.linspace(0, 2 * 1.4848, 32),  # ammonium
}
ablation_result_array = wt_ec.ablation_grid(exch_rate_dict)
# Dump data
with open("ec_grid_glc_amm.pkl", "wb") as handle:
    pickle.dump(ablation_result_array, handle, protocol=pickle.HIGHEST_PROTOCOL)

# ### Pyruvate × ammonium
exch_rate_dict = {
    "r_2033": np.linspace(0, 2 * 4.4444, 32),  # pyruvate
    "r_1654": np.linspace(0, 2 * 1.0, 32),  # ammonium
}
ablation_result_array = wt_ec.ablation_grid(exch_rate_dict)
# Dump data
with open("ec_grid_pyr_amm.pkl", "wb") as handle:
    pickle.dump(ablation_result_array, handle, protocol=pickle.HIGHEST_PROTOCOL)
