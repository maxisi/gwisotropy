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
import pandas as pd

###############################################################################
# LOAD POSTERIORS
###############################################################################

fname = paths.vectors_bbh
with open(fname, 'rb') as f:
    vector_dict = pkl.load(f)
print(f"Loaded: {fname}")

events = sorted(list(vector_dict['n'].keys()))

# pandas gymnastics to get a table with 20 rows and 3 columns,
# filling missing entries with NaN
s = pd.Series(events)
s.index = pd.MultiIndex.from_tuples(s.index.map(lambda x: (x//20, x % 20)))
s = s.unstack(0)

# turn pd.Series into LaTex table
tex = s.to_latex(index=False, header=False, column_format='rrr', na_rep='',
                 caption="Events considered.", label="tab:events")

with open(paths.event_table, 'w') as f:
    f.write(tex)
print(f"Saved: {paths.event_table}")
