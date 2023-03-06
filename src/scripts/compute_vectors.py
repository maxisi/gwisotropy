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
from tqdm import tqdm
from parse import parse
import numpy as np
import utils
import utils.pops 
import utils.transf
import pickle as pkl
import h5py

RNG = np.random.default_rng(12345)

###############################################################################
# LOAD POSTERIORS
###############################################################################

pe_path = str(paths.pe / 'IGWN-{version}-GW{eid}_PEDataRelease_mixed_cosmo.h5')
pe_paths = glob(pe_path.format(version='*', eid='*'))

samples_dict = {}
for path in tqdm(sorted(pe_paths)):
    event = 'GW' + parse(pe_path, path)['eid']
    # select events detected in O3 (i.e., in 2019 or 2020)
    if 'GW19' in event or 'GW20' in event:
        try:
            with h5py.File(path, 'r') as f:
                if 'PublicationSamples' in f:
                    # choose preferred samples
                    k = 'PublicationSamples/posterior_samples'
                    s = f[k][()]
                else:
                    # look for production IMRPhenomXPH run
                    ks = [k for k in f.keys() if ('Prod' in k) or 
                          (('PhenomXPHM' in k) and not ('comoving' in k))]
                    if len(ks) == 0:
                        # no matching runs found
                        print(f.keys())
                        k = None
                        break
                    else:
                        k = ks[0]
                        s = f[k]['posterior_samples'][()]
                samples_dict[event] = s
        except OSError:
            print(f"Unable to load: {path}")

# -----------------------------------------------------------------------------
# throw out events that might include a component below the minimum mass we are
# allowing (M_MIN)

samples_dict_all = samples_dict.copy()
samples_dict = {e: s for e,s in samples_dict_all.items() if
                np.quantile(s['mass_2_source'], 0.01) > utils.MMIN}


###############################################################################
# REWEIGHT POSTERIORS AND FIX POLARIZATION
###############################################################################

# The samples we loaded assumed a prior unfirm in comoving volume and uniform
# in $(1+z)m_1$ (in 2D), we will have to reweight to the MD star formation rate
# and the desired mass population.

Nsamp = utils.NSAMP

k1s = ['mass_1_source', 'mass_ratio', 'redshift']
k2s = k1s + ['ra', 'dec', 'iota', 'psi']

rwt_samples_dict = {}
for e, x in tqdm(samples_dict.items()):
    wts = utils.pops.mass_redshift_pop_wt(*[x[k] for k in k1s]) /\
          utils.pops.li_prior_wt(*[x[k] for k in k2s], cosmo_prior=True)
    rwt_samples_dict[e] = RNG.choice(x, p=wts/sum(wts), size=Nsamp)


# -----------------------------------------------------------------------------
# fix polarization angle

# Before proceeding, we need to fix the $\psi$ samples. This parameter only
# affects the waveform as $2 \psi$, so LALInference only samples it the range
# $\psi \in [0, \pi]$---the posterior for $\psi \in [\pi, 2\pi]$ should be
# identical. RIFT, however, allows for the full range.

rwt_samples_dict = {e: utils.transf.fix_psi(s) for e,s in rwt_samples_dict.items()}


###############################################################################
# GET LOCATION AND ORIENTATION VECTORS
###############################################################################

bbh_vecs_n, bbh_vecs_l, bbh_vecs_j = utils.transf.get_vectors(rwt_samples_dict)
vector_dict = {k: dict(zip(rwt_samples_dict.keys(), vs)) for k, vs in
               zip('nlj', [bbh_vecs_n, bbh_vecs_l, bbh_vecs_j])}

fname = paths.vectors_bbh
with open(fname, 'wb') as f:
    pkl.dump(vector_dict, f, protocol=-1)
print(f"Saved: {fname}")

