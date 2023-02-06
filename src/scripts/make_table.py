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
