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
# LOAD FIT RESULTS
###############################################################################

# number of event number doublings
niter = 5

fit_path = str(paths.data / "control_isotropized/fit_{}.nc")
incr_fits = [az.from_netcdf(fit_path.format(i)) for i in range(niter)]

nsamp = 1000
keys = [r'$\hat{J}_x$', r'$\hat{J}_y$', r'$\hat{J}_z$']

df = pd.DataFrame()
for i, fake_fit in enumerate(incr_fits):
    x = fake_fit.posterior.vL.values
    df_loc = pd.DataFrame(x.reshape(np.prod(x.shape[:2]), 3), columns=keys)
    df_loc['N'] = i
    df = pd.concat([df, df_loc.sample(nsamp, random_state=rng)],
                   ignore_index=True)

###############################################################################
# PLOT
###############################################################################

g = sns.PairGrid(df, hue='N', diag_sharey=False, corner=True,
                 palette='crest_r')
g.map_lower(sns.kdeplot, levels=[1-0.9], alpha=0.7)
g.map_diag(sns.kdeplot, alpha=0.7)

for i, axs in enumerate(g.axes):
    for j, ax in enumerate(axs):
        if ax is not None:
            ax.axvline(0, ls=':', c='k')
            if not i == j:
                ax.axhline(0, ls=':', c='k')

fname = paths.figures / "control_isotropized.pdf"
g.fig.savefig(fname, bbox_inches="tight")
