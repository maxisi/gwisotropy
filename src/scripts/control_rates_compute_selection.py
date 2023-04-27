# /usr/bin/env python

# -*- coding: utf-8 -*-
#
#       Copyright 2023
#       Maximiliano Isi <max.isi@ligo.org>
#       Will M. Farr <will.farr@ligo.org>
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

import paths
from glob import glob
import os
import tqdm
import numpy as np
import utils
import utils.pops 
import utils.transf
import pandas as pd
import h5py

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
# we don't need to worry about the injected spins because the distribution
# was uniform in a1/a2 and isotropic (so falt in cos_tilt)

df_sel['a1'] = np.sqrt(df_sel['spin_1x']**2 + df_sel['spin_1y']**2 + df_sel['spin_1z']**2)
df_sel['a2'] = np.sqrt(df_sel['spin_2x']**2 + df_sel['spin_2y']**2 + df_sel['spin_2z']**2)

df_sel['cos_tilt_1'] = df_sel['spin_1z'] / df_sel['a1']
df_sel['cos_tilt_2'] = df_sel['spin_2z'] / df_sel['a2']

rpdf = pd.read_csv(paths.rates_ref, sep=' ', header=None,
                   index_col=0).squeeze('columns')
our_wt = utils.pops.rp_wt(df_sel['m1'], df_sel['q'], df_sel['z'],
                          df_sel['a1'], df_sel['a2'],
                          df_sel['cos_tilt_1'], df_sel['cos_tilt_2'],
                          refdf=rpdf)

pdraw = sel_wt_total / our_wt
pdraw[~np.isfinite(pdraw)] = np.inf
df_sel['pdrawangle'] = pdraw

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

fname = paths.rates_vectors_sel
with pd.HDFStore(fname) as store:
    store.put("selection", df_sel_vec)
    store.get_storer("selection").attrs.Ndraw = Ndraw
print(f"Saved: {fname}")

