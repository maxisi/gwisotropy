#! /usr/bin/env python

import paths
import pickle as pkl
import arviz as az
import numpy as np
from utils.prior import draw_prior

RNG = np.random.default_rng(12345)

macros = []

###############################################################################
# LOAD POSTERIORS
###############################################################################

fname = paths.vectors_bbh
with open(fname, 'rb') as f:
    vector_dict = pkl.load(f)
print(f"Loaded: {fname}")

Nevents = len(vector_dict['n'])

macros.append("\\renewcommand{\\Nevents}{%i\\xspace}" % Nevents)


###############################################################################
# LOAD FIT
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

import utils.inference as ui

#vdict['prior'] = draw_prior(ndraw=100000, rng=RNG)
for k, vs in vdict.items():
    cl = np.ceil(ui.cl_origin(vs)*100)
    macros.append("\\renewcommand{\\cl%s}{%.0f\\%%\\xspace}" % (k.title(), cl))

###############################################################################
# SAVE MACROS
###############################################################################

with open(paths.macros, 'w') as f:
    f.write("\n".join(macros))
print(f"Saved: {paths.macros}")

