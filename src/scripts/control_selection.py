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
from tqdm import tqdm
import numpy as np
import utils.transf as ut
import utils.inference as ui
import pandas as pd
import pickle as pkl
import arviz as az
import h5py
import pymc as pm
import os

RNG = np.random.default_rng(12345)

###############################################################################
# LOAD PE AND SELECTION VECTORS
###############################################################################

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

fit_path = str(paths.data / "control_selection/fit_{}.nc")
fit_dir = os.path.dirname(fit_path)
if not os.path.exists(fit_dir):
    os.makedirs(fit_dir)

# number of event number doublings
niter = 10

print(f"Running {niter} hierarchical fits---this might take a while!")

rolling_fake_vecs_n_stack = RNG.choice(sel_vecs_n_stack).reshape(1,1,3)
rolling_fake_vecs_j_stack = RNG.choice(sel_vecs_j_stack).reshape(1,1,3)

for i in tqdm(range(niter)):
    model = ui.make_model(rolling_fake_vecs_n_stack, rolling_fake_vecs_j_stack,
                          sel_vecs_n_stack, sel_vecs_j_stack,
                          df_sel_vec.pdrawangle.values, Ndraw)

    with model:
        f = az.convert_to_inference_data(pm.sample())

    fname = fit_path.format(i)
    f.to_netcdf(fname)
    print(f"Saved: {fname}")

    # update running list of fake vectors, doubling the number of entries
    N = fake_vecs_n_stack.shape[0]
    fake_vecs_n_stack = RNG.choice(sel_vecs_n_stack, N).reshape(N,1,3)
    fake_vecs_j_stack = RNG.choice(sel_vecs_j_stack, N).reshape(N,1,3)

    rolling_fake_vecs_n_stack = np.concatenate([rolling_fake_vecs_n_stack,
                                                fake_vecs_n_stack])
    rolling_fake_vecs_j_stack = np.concatenate([rolling_fake_vecs_j_stack,
                                                fake_vecs_j_stack])

