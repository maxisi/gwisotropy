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

def random_shift(vec_stack, rng=RNG):
    """Randomize events over the sky by rotating posteriors rigidly.
    """
    new_vec_stack = []
    # iterate over single-event posteriors
    for vecs in vec_stack:

        # choose random direction to rotate about
        theta = np.arccos(rng.uniform(-1, 1))
        phi = rng.uniform(0, 2*np.pi)

        # choose a random angle to rotate by
        angle = rng.uniform(0, 2*np.pi)

        # create corresponding rotation matrix
        direction = np.array([np.sin(theta)*np.cos(phi),
                              np.sin(theta)*np.sin(phi),
                              np.cos(theta)]).T

        R = ut.rotation_matrix(angle, direction)
        new_vec_stack.append(np.dot(R, vecs.T).T)
    return new_vec_stack

fit_path = str(paths.data / "control_isotropized/fit_{}.nc")
fit_dir = os.path.dirname(fit_path)
if not os.path.exists(fit_dir):
    os.makedirs(fit_dir)

# number of event number doublings
niter = 6

print(f"Running {niter} hierarchical fits---this might take a while!")

for i in tqdm(range(niter)):
    if i == 0:
        rolling_fake_vecs_n_stack = random_shift(bbh_vecs_n_stack)
        rolling_fake_vecs_j_stack = random_shift(bbh_vecs_j_stack)
    else:
        # create new stack of same length as current stack
        fake_vecs_n_stack = random_shift(rolling_fake_vecs_n_stack)
        fake_vecs_j_stack = random_shift(rolling_fake_vecs_j_stack)

        # append to running stack, doubling the number of events
        rolling_fake_vecs_n_stack = np.concatenate([rolling_fake_vecs_n_stack,
                                                    fake_vecs_n_stack])
        rolling_fake_vecs_j_stack = np.concatenate([rolling_fake_vecs_j_stack,
                                                    fake_vecs_j_stack])

    model = ui.make_model(rolling_fake_vecs_n_stack, rolling_fake_vecs_j_stack,
                          sel_vecs_n_stack, sel_vecs_j_stack,
                          df_sel_vec.pdrawangle.values, Ndraw)

    with model:
        f = az.convert_to_inference_data(pm.sample())
    fname = fit_path.format(i)
    f.to_netcdf(fname)
    print(f"Saved: {fname}")
