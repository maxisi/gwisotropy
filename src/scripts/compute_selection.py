import paths
from glob import glob
import os
import tqdm
import numpy as np
import utils
import utils.pops 
import pandas as pd

RNG = np.random.default_rng(12345)

###############################################################################
# LOAD SELECTION DATA
###############################################################################

# load results of O3 injection campaign made public by LIGO-Virgo

with h5py.File(paths.selection, 'r') as f:
    d = {k: np.array(f['injections'][k]) for k in f['injections'].keys()}
    o3_inj_data = pd.DataFrame(d)
    Ndraw = f.attrs['total_generated']

# -----------------------------------------------------------------------------
# rename some variables to work with the code below

rename_dict = {
    'mass1_source': 'm1',
    'redshift': 'z',
    'right_ascension': 'ra',
    'polarization': 'psi'
}
rename_dict.update({f'spin{i}{x}':f'spin_{i}{x}' for i in [1,2] for x in 'xyz'})
df_sel = o3_inj_data.rename(columns=rename_dict)

df_sel['sindec'] = np.sin(df_sel['declination'])
df_sel['q'] = df_sel['mass2_source']/df_sel['m1']
df_sel['cosiota'] = np.cos(df_sel['inclination'])

# -----------------------------------------------------------------------------
# identify detections, defined as a FAR < FARMAX in any of the pipelines

far_keys = ['far_cwb', 'far_gstlal', 'far_mbta', 'far_pycbc_bbh',
            'far_pycbc_hyperbank']
detected_mask = False
for k in far_keys:
    detected_mask = (o3_inj_data[k] < utils.FARMAX_S) | detected_mask

print(f"There were {sum(detected_mask)} detected injections in the set.")


###############################################################################
# COMPUTE WEIGHTS
###############################################################################

# -----------------------------------------------------------------------------
# compute selection weights based on the population sampled for injections

sel_wt_orbit = df_sel['polarization_sampling_pdf'] *\
               df_sel['inclination_sampling_pdf']/np.sin(df_sel['inclination'])

sel_wt_sky = df_sel['declination_sampling_pdf']/np.cos(df_sel['declination'])*\
             df_sel['right_ascension_sampling_pdf']

sel_wt_m1q = df_sel['mass1_source_mass2_source_sampling_pdf']*df_sel['m1']

sel_wt_redshift = df_sel['redshift_sampling_pdf']

sel_wt_total = sel_wt_orbit*sel_wt_sky * sel_wt_m1q * sel_wt_redshift

# -----------------------------------------------------------------------------
# readjust weight to our chosen mass and redshift population

our_wt = utils.pops.mass_redshift_pop_wt(df_sel['m1'], df_sel['q'], df_sel['z'])

df_sel['pdrawangle'] = sel_wt_total / our_wt

df_sel_raw = df_sel.copy()
df_sel = df_sel_raw[detected_mask].copy()


###############################################################################
# GET LOCATION AND ORIENTATION VECTORS
###############################################################################

sel_vecs_n, _, sel_vecs_j = utils.transf.get_sel_vectors(df_sel)

# form new DataFrame with vectors and selection
df_n = pd.DataFrame([xyz for xyz in sel_vecs_n.values],
                    columns=['nx', 'ny', 'nz'], index=sel_vecs_n.index)
df_j = pd.DataFrame([xyz for xyz in sel_vecs_j.values],
                    columns=['jx','jy', 'jz'], index=sel_vecs_j.index)

df_sel_vec = pd.concat([df_n, df_j], axis=1)
df_sel_vec['pdrawangle'] = df_sel['pdrawangle']

fname = paths.vectors_sel
with pd.HDFStore(fname) as store:
    store.put("selection", df_sel_vec)
    store.get_storer("selection").attrs.Ndraw = Ndraw
print(f"Saved: {fname}")

# vector_dict = dict(zip('nlj', [sel_vecs_n, sel_vecs_l, sel_vecs_j]))
# 
# fname = paths.vectors_sel
# with open(fname, 'wb') as f:
#     pkl.dump(vector_dict, f, protocol=-1)
# print(f"Saved: {fname}")

