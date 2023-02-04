import paths
import tqdm
import numpy as np
import utils
import utils.inference as ui
import pandas as pd
import pickle as pkl
import arviz as az
import h5py
import pymc as pm

RNG = np.random.default_rng(12345)

###############################################################################
# LOAD PE AND SELECTION VECTORS
###############################################################################

# -----------------------------------------------------------------------------
# retrieve location and orientation vectors for observed BBHs

with open(paths.vectors_bbh, 'rb') as f:
    vector_dict_bbh = pkl.load(f)

bbh_vecs_n_stack = np.stack(list(vector_dict_bbh['n'].values()))
bbh_vecs_j_stack = np.stack(list(vector_dict_bbh['j'].values()))

# -----------------------------------------------------------------------------
# retrieve location and orientation vectors for detected injections

df_sel_vec = pd.read_hdf(paths.vectors_sel, "selection")
with h5py.File(paths.vectors_sel) as f:
    Ndraw = f['selection'].attrs['Ndraw']

sel_vecs_n_stack = df_sel_vec[[f'n{i}' for i in 'xyz']].values
sel_vecs_j_stack = df_sel_vec[[f'j{i}' for i in 'xyz']].values

###############################################################################
# SAMPLE
###############################################################################

model = ui.make_model(bbh_vecs_n_stack, bbh_vecs_j_stack,
                      sel_vecs_n_stack, sel_vecs_j_stack,
                      df_sel_vec.pdrawangle.values, Ndraw)
with model:
    trace = pm.sample()
    fit = az.convert_to_inference_data(trace)

fname = paths.result
fit.to_netcdf(fname)
print(f"Saved: {fname}")
