# /usr/bin/env python

"""Sets up a slurm workflow to produce individual-event skymaps included in
static/ directory. This requires the ligo.skymap pakcage
(https://lscsoft.docs.ligo.org/ligo.skymap/index.html), which is not a
otherwise a dependency for this paper.
"""

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
import pickle as pkl
from astropy.table import Table
import ligo.skymap.io
import os

###############################################################################
# LOAD PE AND SELECTION VECTORS
###############################################################################

# -----------------------------------------------------------------------------
# retrieve location and orientation vectors for observed BBHs

with open(paths.vectors_bbh, 'rb') as f:
    vector_dict_bbh = pkl.load(f)

bbh_vecs_j = list(vector_dict_bbh['j'].values())
events = list(vector_dict_bbh['j'].keys())

###############################################################################
# PROCESS
###############################################################################

# export J-vector coordinate samples for all events

def export_data(keys, vectors, path=None, **kwargs):
    idxs = ...
    for key, vecs in zip(keys, tqdm(vectors)):
        dist = np.linalg.norm(vecs, axis=1)
        ra = np.arctan2(vecs[idxs,1], vecs[idxs,0])
        dec = np.arcsin(vecs[idxs,2]/dist)
        table = Table([ra, dec, dist], names=['ra', 'dec', 'dist'])
        p = (path or "{key}.hdf5").format(key=key)
        bdir = os.path.dirname(p)
        if not os.path.exists(bdir):
            os.makedirs(bdir)
        ligo.skymap.io.write_samples(table, p, path='posterior_samples', **kwargs)

outpath = str(paths.data / 'cache/posterior_j_{key}.h5')
os.makedirs(os.path.dirname(outpath), exist_ok=True)
export_data(events, bbh_vecs_j, outpath, overwrite=True)


# create a slurm disBatch TaskFile to obtain a fits file from the above samples
fitsdir = str(paths.data / 'fits')
os.makedirs(fitsdir, exist_ok=True)

logpath = str(paths.data / 'fits/logs/{key}.log')
os.makedirs(os.path.dirname(logpath), exist_ok=True)

l = "ligo-skymap-from-samples -j 20 -o %s --fitsoutname skymap_j_{key}.fits.gz %s &> %s" % (fitsdir, outpath, logpath)

lines = [l.format(key=k) for k in events]
with open(paths.data / "fits/FitsTaskFile", 'w') as f:
    f.write('\n'.join(lines))

# create a slurm disBatch TaskFile to plot skymaps
figdir = paths.static
os.makedirs(figdir, exist_ok=True)

l = "ligo-skymap-plot %s/skymap_j_{key}.fits.gz -o %s/skymap_j_{key}.{ext} --figure-height 4 &> %s_{ext}.log" % (fitsdir, figdir, logpath)

lines = [l.format(key=k, ext='pdf') for k in events]
with open(paths.data / "fits/PlotFitsTaskFile", 'w') as f:
    f.write('\n'.join(lines))
