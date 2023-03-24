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
import pickle as pkl
import arviz as az
import numpy as np
import utils
from utils.prior import draw_prior
import utils.inference as ui

RNG = np.random.default_rng(12345)

macros = []

###############################################################################
# INDIVIDUAL EVENTS
###############################################################################

fname = paths.vectors_bbh
with open(fname, 'rb') as f:
    vector_dict = pkl.load(f)
print(f"Loaded: {fname}")

Nevents = len(vector_dict['n'])

macros.append("\\renewcommand{\\Nevents}{%i\\xspace}" % Nevents)


###############################################################################
# HIERARCHICAL FIT
###############################################################################

fit = az.from_netcdf(paths.result)

x = fit.posterior.vN.values
nhats = x.reshape(np.prod(x.shape[:2]), 3)

x = fit.posterior.vL.values
jhats = x.reshape(np.prod(x.shape[:2]), 3)

# vprior = draw_prior(ndraw=100000, rng=RNG)
std_prior = 0.3089

vdict = {'n': nhats, 'j': jhats}

for k, vs in vdict.items():
    std = np.std(vs, axis=0)
    x = (std - std_prior)/std_prior
    for i, xi in zip('xyz', x):
        macros.append("\\renewcommand{\\varimp%s%s}{%.0f\\%%\\xspace}"
                      % (k.upper(),i,np.abs(xi)*100))


#vdict['prior'] = draw_prior(ndraw=100000, rng=RNG)
for k, vs in vdict.items():
    cl = np.ceil(ui.cl_origin(vs)*100)
    macros.append("\\renewcommand{\\cl%s}{%.0f\\%%\\xspace}" % (k.title(), cl))

###############################################################################
# VALIDATION RUNS
###############################################################################

macros.append(r"\renewcommand{\Niteriso}{%i\xspace}" % (utils.NITER_ISO-1))
nmaxiso = Nevents * 2**(utils.NITER_ISO - 1)
macros.append(r"\renewcommand{\Nmaxiso}{%i\xspace}" % nmaxiso)

macros.append(r"\renewcommand{\Nitersel}{%i\xspace}" % (utils.NITER_SEL-1))
macros.append(r"\renewcommand{\Nstartsel}{%i\xspace}" % utils.NSTART_SEL)
nmaxsel = utils.NSTART_SEL * 2**(utils.NITER_SEL - 1)
macros.append(r"\renewcommand{\Nmaxsel}{%i\xspace}" % nmaxsel)

###############################################################################
# SAVE MACROS
###############################################################################

with open(paths.macros, 'w') as f:
    f.write("\n".join(macros))
print(f"Saved: {paths.macros}")

