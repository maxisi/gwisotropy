"""
Exposes common paths useful for manipulating datasets and generating figures.

"""
from pathlib import Path

# Absolute path to the top level of the repository
root = Path(__file__).resolve().parents[2].absolute()

# Absolute path to the `src` folder
src = root / "src"

# Absolute path to the `src/data` folder (contains datasets)
data = src / "data"
pe = data / "pe_gwosc_o3"
vectors_bbh = data / "vectors_bbh.pkl"
selection = data / "endo3_bbhpop-LIGO-T2100113-v12.hdf5"
vectors_sel = data / "vectors_sel.hdf5"
result = data / "gwisotropy_result.nc"

# Absolute path to the `src/static` folder (contains static images)
static = src / "static"

# Absolute path to the `src/scripts` folder (contains figure/pipeline scripts)
scripts = src / "scripts"

# Absolute path to the `src/tex` folder (contains the manuscript)
tex = src / "tex"

# Absolute path to the `src/tex/figures` folder (contains figure output)
figures = tex / "figures"

# Absolute path to the `src/tex/output` folder (contains other user-defined output)
output = tex / "output"
