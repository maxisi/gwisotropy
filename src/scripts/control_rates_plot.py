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
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import arviz as az
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.lines import Line2D

sns.set(context='notebook', palette='colorblind', font_scale=2)

RNG = np.random.default_rng(12345)


###############################################################################
# LOAD FIT
###############################################################################

# main result
fit = az.from_netcdf(paths.result)

x = fit.posterior.vN.values
nhats_main = x.reshape(np.prod(x.shape[:2]), 3)

x = fit.posterior.vL.values
jhats_main = x.reshape(np.prod(x.shape[:2]), 3)

# control result
fit = az.from_netcdf(paths.rates_result)

x = fit.posterior.vN.values
nhats_rates = x.reshape(np.prod(x.shape[:2]), 3)

x = fit.posterior.vL.values
jhats_rates = x.reshape(np.prod(x.shape[:2]), 3)

###############################################################################
# PLOT
###############################################################################

ckeys = [r'$v_{N,x}$', r'$v_{N,y}$', r'$v_{N,z}$',
         r'$v_{J,x}$', r'$v_{J,y}$', r'$v_{J,z}$']
df0 = pd.DataFrame(np.hstack((nhats_main, jhats_main)), columns=ckeys)
df1 = pd.DataFrame(np.hstack((nhats_rates, jhats_rates)), columns=ckeys)

c1 = 'k'
lkws = dict(ls=':', c='k')
with sns.axes_style("ticks"):
    pg = sns.PairGrid(df0, corner=True)#, hue='run', hue_kws={'fill': [True, False]})
    pg.map_diag(sns.kdeplot, fill=True, common_norm=False)
    pg.map_lower(sns.kdeplot, fill=True, alpha=0.7, thresh=0.1)

    pg.data = df1
    pg.map_diag(sns.kdeplot, fill=False, color=c1, common_norm=False)
    pg.map_lower(sns.kdeplot, fill=False, color=c1, levels=[0.1, 0.5, 0.9],
                 thresh=0.1)

    for i, axs in enumerate(pg.axes):
        for j, ax in enumerate(axs):
            if ax:
                ax.axvline(0, **lkws)
                ax.set_xlim(-1, 1)
                if i != j:
                    ax.axhline(0, **lkws)
                    ax.set_ylim(-1, 1)
                    if j == 0:
                        ax.set_yticks([-1, 0, 1])

p = paths.figures / "control_rates_jn_corner.pdf"
pg.fig.savefig(p, bbox_inches="tight", dpi=300)
print(f"Saved: {p}")
