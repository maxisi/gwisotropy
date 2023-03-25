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
import utils.plots
import pandas as pd
import arviz as az
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.lines import Line2D
from matplotlib.figure import Figure
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

sns.set(context='notebook', palette='colorblind', font_scale=1.5)

RNG = np.random.default_rng(12345)


###############################################################################
# LOAD FIT
###############################################################################

fit = az.from_netcdf(paths.result)

x = fit.posterior.vN.values
nhats = x.reshape(np.prod(x.shape[:2]), 3)

x = fit.posterior.vL.values
jhats = x.reshape(np.prod(x.shape[:2]), 3)

###############################################################################
# PLOT
###############################################################################

ckeys = [r'$v_{N,x}$', r'$v_{N,y}$', r'$v_{N,z}$',
         r'$v_{J,x}$', r'$v_{J,y}$', r'$v_{J,z}$']
df = pd.DataFrame(np.hstack((nhats, jhats)), columns=ckeys)

lkws = dict(ls=':', c='k')
with sns.axes_style("ticks"):
    pg = sns.PairGrid(df, corner=True)
    pg.map_diag(sns.kdeplot, fill=True)
    pg.map_lower(sns.kdeplot, fill=True, thresh=0.1)
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

akws = dict(arrowstyle="->", color=sns.color_palette()[0], alpha=0.3, lw=2.5)

axins = pg.axes[1,1].inset_axes([3.45, 0.5, 4, 2])
arr_img = plt.imread(paths.static / "dipole_skymap_n.png")
axins.imshow(arr_img)
axins.set_title(r"$\hat{v}_N$", pad=-20)
axins.axis('off')
axins.annotate("", xy=(-0.75, 0.5), xytext=(-0.25, 0.5),
               arrowprops=akws, xycoords='axes fraction')

axins = pg.axes[4,4].inset_axes([-0.25, 2, 4, 2])
arr_img = plt.imread(paths.static / "dipole_skymap_j.png")
axins.imshow(arr_img)
axins.set_title(r"$\hat{v}_J$", pad=-20)
axins.axis('off')
axins.annotate("", xy=(0.5, -0.75), xytext=(0.5, -0.25),
           arrowprops=akws, xycoords='axes fraction')

xys = [([-0.02, -0.02], [0.52, 0.95]), ([0.52, 1], [0, 0])]
lines = [Line2D(x, y, lw=5., alpha=0.3, figure=pg.fig, transform=pg.fig.transFigure,
               ) for x, y in xys]
pg.fig.patches.extend(lines)

p = paths.figures / "jn_corner.pdf"
pg.fig.savefig(p, bbox_inches="tight", dpi=300)
print(f"Saved: {p}")
