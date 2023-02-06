#! /usr/bin/env python

import paths
import pickle as pkl

###############################################################################
# LOAD POSTERIORS
###############################################################################

fname = paths.vectors_bbh
with open(fname, 'rb') as f:
    vector_dict = pkl.load(f)
print(f"Loaded: {fname}")

Nevents = len(vector_dict['n'])

macros = ["\\renewcommand{\\Nevents}{%i\\xspace}" % Nevents]

with open(paths.macros, 'w') as f:
    f.write("\n".join(macros))
print(f"Saved: {paths.macros}")
