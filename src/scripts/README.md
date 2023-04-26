# Software

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
    classDef data fill:#fbf
    S --> F((hierarchical_fit.py)):::script;
    B[(GWOSC)]--> |download_gwosc_pe.py|PE(pe_gwosc_o3/):::data;
    PE --> compute_vectors.py;
    compute_vectors.py --> D(vectors_bbh.pkl):::data;
    D --> F;
    F --> R(gwisotropy_result.nc):::data;
    F2a{{dipole_skymap_j.png}}:::fig --> corner_plot.py;
    classDef fig fill:#f96
    F2b{{dipole_skymap_n.png}}:::fig --> corner_plot.py;
    SKY[/setup_dipole_skymaps.py/] -.-> F2a;
    SKY -.-> F2b;
    R --> corner_plot.py;
    corner_plot.py --> F2{{jn_corner.pdf}}:::fig;
    F2 --> M{manuscript};
    R --> density_3d_plot.py;
    density_3d_plot.py --> F3a{{density_3d_j.pdf}}:::fig;
    density_3d_plot.py --> F3b{{density_3d_n.pdf}}:::fig;
    density_3d_plot.py --> F3c{{density_3d_prior.pdf}}:::fig;
    R -.-> SKY;
    F3a --> M;
    F3b --> M;
    F3c --> M;
    R --> norm_plot.py;
    norm_plot.py --> F4{{jn_norm.pdf}}:::fig;
    F4 --> M;
    R --> make_macros.py;
    make_macros.py --> m>macros.tex];
    m --> M;
    R --> make_table.py;
    make_table.py --> e>events.tex];
    e --> M;
    NB[/gwiso_vector_plots.nb/] -.-> F1a{{gw_orientation_iso.pdf}}:::fig;
    NB -.-> F1b{{gw_orientation_dip.pdf}}:::fig;
    NB -.-> F1c{{gw_orientation_inc.pdf}}:::fig;
    F1a --> M;
    F1b --> M;
    F1c --> M;
    prior_plot.py --> P{{prior.pdf}}:::fig;
    P --> M;
    D --> control_selection.py;
    S --> control_selection.py;
    control_selection.py --> CSEL(control_selection/):::data;
    CSEL --> plot_control_selection.py;
    plot_control_selection.py --> CSELPLOT{{control_selection.pdf}}:::fig;
    CSELPLOT --> M;
    SKYALL[/setup_event_skymaps.py/] -.-> SMAPS{{skymap_j_*.pdf}}:::fig
    D -.-> SKYALL;
    SMAPS --> s>skymaps.tex];
    make_skymap_figure.py --> s;
    s --> M;
    X[(10.5281/zenodo.7843926)]-->|control_rates_get_lvk_pop.py|RP(control_rates_powerlawpeak_map.txt):::data;
    RP --> control_rates_compute_selection.py;
    control_rates_compute_selection.py --> CS(control_rates_vectors_sel.hdf5):::data;
    CS --> CF(control_rates_hierarchical_fit.py):::script
    RP --> control_rates_compute_vectors.py;
    control_rates_compute_vectors.py --> CD(control_rates_vectors_bbh.pkl):::data;
    CD --> CF;
    CF --> CR(control_rates_gwisotropy_result.nc):::data;
    CR --> control_rates_plot.py;
    R --> control_rates_plot.py;
    control_rates_plot.py --> CF2{{contro_rates_jn_corner.pdf}}:::fig;
    CF2 --> M{manuscript};
```

The scripts rely on a small Python package provided in the [`utils`](utils) directory, which contains its own documentation.

The environment requirements to execute these scripts are specified in the [`environment.yml`](https://github.com/maxisi/gwisotropy/blob/main/environment.yml) in the source directory; if not using `showyourwork`, this can be used to create a compatible Conda environment by doing, e.g.,

```
conda env create -f environment.yml
```

Some intermediate data products are generated in the process and cached by
_showyourwork_ to speed up computation and skip some intermediate steps. The cached data products will be automatically downloaded from [10.5281/zenodo.7775266](https://doi.org/10.5281/zenodo.7775266) when using _showyourwork_.
