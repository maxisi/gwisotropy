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

rule table:
    input:
        "src/data/vectors_bbh.pkl"
    output:
        "src/tex/output/events.tex"
    script:
        "src/scripts/make_table.py"

rule macros:
    input:
        "src/data/vectors_bbh.pkl"
    output:
        "src/tex/output/macros.tex"
    script:
        "src/scripts/make_macros.py"
