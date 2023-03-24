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
import pandas as pd
import pickle as pkl
import arviz as az
import h5py
import pymc as pm
import seaborn as sns
from glob import glob
import matplotlib.pyplot as plt
import utils.settings

sns.set(context='notebook', palette='colorblind', font_scale=1.5)
RNG = np.random.default_rng(12345)

###############################################################################
# LOAD FIT RESULTS
###############################################################################

fit_path = str(paths.data / "control_selection/fit_{}.nc")
nfits = len(glob(fit_path.format('*')))
incr_fits = [az.from_netcdf(fit_path.format(i)) for i in range(nfits)]

nsamp = 4000
keys = [r'$\vec{v}_{J,x}$', r'$\vec{v}_{J,y}$', r'$\vec{v}_{J,z}$']

df = pd.DataFrame()
for i, fake_fit in enumerate(incr_fits):
    x = fake_fit.posterior.vL.values
    df_loc = pd.DataFrame(x.reshape(np.prod(x.shape[:2]), 3), columns=keys)
    df_loc['N'] = i
    df = pd.concat([df, df_loc.sample(nsamp, random_state=RNG)],
                   ignore_index=True)

###############################################################################
# PLOT
###############################################################################

kde_kws = dict(common_norm=False, alpha=0.7)
cmap = 'crest_r'
g = sns.PairGrid(df, hue='N', diag_sharey=False, corner=True, palette=cmap)
g.map_lower(sns.kdeplot, levels=1, thresh=0.1, **kde_kws)
g.map_diag(sns.kdeplot, **kde_kws)

for i, axs in enumerate(g.axes):
    for j, ax in enumerate(axs):
        if ax is not None:
            ax.axvline(0, ls=':', c='k')
            ax.set_xlim(-1, 1)
            if not i == j:
                ax.axhline(0, ls=':', c='k')
                ax.set_ylim(-1, 1)
                if j == 0:
                    ax.set_yticks([-1, 0, 1])

# backtrack number of detections per iteration
nstart = utils.settings.NSTART_SEL
niter = utils.settings.NITER_SEL
ns = [nstart*2**i for i in range(niter)]

# add colorbar
norm = plt.Normalize(vmin=min(ns), vmax=max(ns))
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

cax = g.fig.add_axes([g.axes[1,1].get_position().x0,
                      g.axes[0,0].get_position().y0,
                      g.axes[0,0].get_position().width*0.1,
                      g.axes[0,0].get_position().height])
plt.colorbar(sm, cax=cax, orientation='vertical')
yt = cax.get_yticks()
cax.set_yticklabels([f'$N = {n:.0f}$' for n in yt], fontsize=14)

# add norm inset
ax0 = g.axes[0,0]
ax1 = g.axes[-1,-1]
ax = g.fig.add_axes([ax1.get_position().x0, ax0.get_position().y0,
                     ax0.get_position().width, ax0.get_position().height])

df['norm'] = np.linalg.norm([df[k] for k in keys], axis=0)
sns.kdeplot(data=df[['norm', 'N']], x='norm', hue='N', ax=ax, palette=cmap,
            common_norm=False)

ax.get_legend().remove()
ax.set_xlabel(r"$|\vec{v}_J|$")
ax.set_ylabel('')
ax.get_yaxis().set_visible(False)
ax.set_xlim(0, 1)


fname = paths.figures / "control_selection.pdf"
g.fig.savefig(fname, bbox_inches="tight")
