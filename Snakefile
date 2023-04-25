rule gwosc:
    output:
        directory("src/data/pe_gwosc_o3")
    script:
        "src/scripts/download_gwosc_pe.py"

rule vectors:
    input:
        "src/data/pe_gwosc_o3"
    output:
        "src/data/vectors_bbh.pkl"
    cache:
        True
    script:
        "src/scripts/compute_vectors.py"

rule selection:
    input:
        "src/data/endo3_bbhpop-LIGO-T2100113-v12.hdf5"
    output:
        "src/data/vectors_sel.hdf5"
    cache:
        True
    script:
        "src/scripts/compute_selection.py"

rule fit:
    input:
        "src/data/vectors_sel.hdf5",
        "src/data/vectors_bbh.pkl"
    output:
        "src/data/gwisotropy_result.nc"
    cache:
        True
    script:
        "src/scripts/hierarchical_fit.py"

rule controlsel:
    input:
        "src/data/vectors_sel.hdf5",
    output:
        directory("src/data/control_selection")
    cache:
        True
    script:
        "src/scripts/control_selection.py"

rule controlselplot:
    input:
        "src/data/control_selection"
    output:
        "src/tex/figures/control_selection.pdf"
    script:
        "src/scripts/control_selection_plot.py"

rule table:
    input:
        "src/data/vectors_bbh.pkl"
    output:
        "src/tex/output/events.tex"
    script:
        "src/scripts/make_table.py"

rule skymaps:
    output:
        "src/tex/output/skymaps.tex"
    script:
        "src/scripts/make_skymap_figure.py"

rule macros:
    input:
        "src/data/gwisotropy_result.nc",
        "src/scripts/utils/settings.py",
        "src/data/vectors_bbh.pkl"
    output:
        "src/tex/output/macros.tex"
    script:
        "src/scripts/make_macros.py"

rule contrlratesgetlvk:
    input:
        "src/data/analyses_PowerLawPeak/analyses/PowerLawPeak/o1o2o3_mass_c_iid_mag_iid_tilt_powerlaw_redshift_result.json"
    output:
        "src/data/control_rates_powerlawpeak_map.txt"
    cache:
        True
    script:
        "src/scripts/control_rates_get_lvk_pop.py"

rule controlratesvectors:
    input:
        "src/data/pe_gwosc_o3",
        "src/data/control_rates_powerlawpeak_map.txt"
    output:
        "src/data/control_rates_vectors_bbh.pkl"
    cache:
        True
    script:
        "src/scripts/control_rates_compute_vectors.py"

rule controlratesselection:
    input:
        "src/data/endo3_bbhpop-LIGO-T2100113-v12.hdf5",
        "src/data/control_rates_powerlawpeak_map.txt"
    output:
        "src/data/control_rates_vectors_sel.hdf5"
    cache:
        True
    script:
        "src/scripts/control_rates_compute_selection.py"

rule controlratesfit:
    input:
        "src/data/control_rates_vectors_sel.hdf5",
        "src/data/control_rates_vectors_bbh.pkl"
    output:
        "src/data/control_rates_gwisotropy_result.nc"
    cache:
        True
    script:
        "src/scripts/control_rates_hierarchical_fit.py"

