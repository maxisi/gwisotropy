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
from parse import parse

sns.set(context='notebook', palette='colorblind', font_scale=1.5)
RNG = np.random.default_rng(12345)

###############################################################################
# LOAD FIT RESULTS
###############################################################################

# load fits in order of number of vectors in simulated catalog
fit_path = str(paths.data / "control_isotropized/fit_{}.nc")
nvecs = sorted([int(parse(fit_path,p)[0]) for p in glob(fit_path.format('*'))])
incr_fits = [az.from_netcdf(fit_path.format(n)) for n in nvecs]

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

#with sns.axes_style("whitegrid", {"grid.linestyle": ':'}):
cmap = 'crest_r'
kde_kws = dict(common_norm=False, alpha=0.7)
g = sns.PairGrid(df, hue='N', diag_sharey=False, corner=True, palette=cmap)
g.map_lower(sns.kdeplot, levels=1, thresh=0.1, **kde_kws)
g.map_diag(sns.kdeplot, **kde_kws)

ax = g.axes[0,0]
cs = sns.color_palette(cmap, n_colors=len(nvecs))
for n, c in zip(nvecs, cs):
    ax.plot([], [], label="$N = {}$".format(n), c=c)
ax.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), frameon=False,
          fontsize=14)

for i, axs in enumerate(g.axes):
    for j, ax in enumerate(axs):
        if ax is not None:
            ax.axvline(0, ls=':', c='k')
            ax.set_xlim(-1, 1)
            if i != j:
                ax.axhline(0, ls=':', c='k')
                ax.set_ylim(-1, 1)
                if j == 0:
                    ax.set_yticks([-1, 0, 1])

ax0 = g.axes[0,0]
ax1 = g.axes[-1,-1]
ax = g.fig.add_axes([ax1.get_position().x0, ax0.get_position().y0,
                     ax0.get_position().width, ax0.get_position().height])
df['norm'] = np.linalg.norm([df[k] for k in keys], axis=0)
df_norm = df[['norm', 'N']]
sns.kdeplot(data=df_norm, x='norm', hue='N', ax=ax, palette=cmap, common_norm=False)
ax.get_legend().remove()
ax.set_xlabel(r"$|\vec{v}_J|$")
ax.set_ylabel('')
ax.get_yaxis().set_visible(False)
ax.set_xlim(0, 1)

fname = paths.figures / "control_isotropized.pdf"
g.fig.savefig(fname, bbox_inches="tight")
