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
import seaborn as sns
from utils.prior import draw_prior
import numpy as np

sns.set(context='notebook', palette='colorblind')

RNG = np.random.default_rng(12345)

# draw vectors from prior
v = draw_prior(ndraw=100000, rng=RNG)

# plot 2D slice
with sns.axes_style("ticks"):
    g = sns.jointplot(x=v[:,0], y=v[:,1], kind='hex');
    g.ax_joint.set_xlabel(r"$\vec{v}_x$")
    g.ax_joint.set_ylabel(r"$\vec{v}_y$");

    p = paths.figures / "prior.pdf"
    g.fig.savefig(p, bbox_inches='tight')
    print(f"Saved: {p}")
