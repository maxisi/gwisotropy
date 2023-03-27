# Software

## Plotting scripts

This directory contains the files to reproduce the analysis in _The Directional Isotropy of LIGO-Virgo binaries_ (Isi+2023), 
from data publicly avaliable in Zenodo: parameter estimation data releases from [GWOSC](https://www.gw-openscience.org) and LIGO-Virgo selection function
injections [10.5281/zenodo.5546676](https://doi.org/10.5281/zenodo.5546676); intermediate data products produced in our analysis are cached in
[10.5281/zenodo.7775266](https://doi.org/10.5281/zenodo.7775266).

When compiling the manuscript, _showyourwork_ will automatically use cached data products to speed of the computation.

The relation between the different pieces of software and data is as follows:

```mermaid
graph TD;
    A[(10.5281/zenodo.5546676)]-->compute_selection.py;
    compute_selection.py --> S(vectors_sel.hdf5):::data;
    classDef data fill:#f96
    S --> F[hierarchical_fit.py]:::script;
    B[(GWOSC)]--> download_gwosc_pe.py;
    download_gwosc_pe.py-->PE(pe_gwosc_o3):::data;
    PE --> compute_vectors.py;
    compute_vectors.py --> D(vectors_bbh.pkl):::data;
    D --> F;
    F --> R(gwisotropy_result.nc):::data;
    F2a{{dipole_skymap_j.png}}:::fig --> corner_plot.py;
    F2b{{dipole_skymap_n.png}}:::fig --> corner_plot.py;
    SKY[/setup_dipole_skymaps.py/] -.-> F2a;
    SKY -.-> F2b;
    R --> corner_plot.py;
    corner_plot.py --> F2{{jn_corner.pdf}}:::fig;
    F2 --> M{manuscript};
    R --> density_3d_plot.py;
    density_3d_plot.py --> F3a{{density_3d_j.pdf}};
    density_3d_plot.py --> F3b{{density_3d_n.pdf}};
    density_3d_plot.py --> F3c{{density_3d_prior.pdf}};
    F3a --> M;
    F3b --> M;
    F3c --> M;
    R --> norm_plot.py;
    norm_plot.py --> F4{{jn_norm.pdf}};
    F4 --> M;
```

The environment requirements to execute these scripts are specified in the [`environment.yml`](https://github.com/maxisi/gwisotropy/blob/main/environment.yml) in the source directory; if not using `showyourwork`, this can be used to create a compatible Conda environment by doing, e.g.,

```
conda env create -f environment.yml
```

Some intermediate data products are generated in the process and cached by
_showyourwork_ to speed up computation and skip some intermediate steps.
