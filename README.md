#  The directional isotropy of LIGO-Virgo binaries

Software and data release for _The directional isotropy of LIGO-Virgo binaries_ by M. Isi, W. M. Farr and V. Varma (2023).

👉 The manuscript is programatically generated using _[showyourwork](https://show-your.work/en/latest/)_; scripts and auxiliary code to produce the results and figures are contained in **[`src/scripts`](src/scripts)**, please see that directory for documentation. Data are stored in Zenodo (https://doi.org/10.5281/zenodo.7775266).

## Abstract
We demonstrate how to constrain the degree of absolute alignment of the total angular momenta
of LIGO-Virgo binary black holes, looking for a special direction in space that would break isotropy.
We also allow for inhomogeneities in the distribution of black holes over the sky. Making use of
dipolar models for the spatial distribution and orientation of the sources, we analyze 57 signals with
false-alarm rates ≤ 1/yr from the third LIGO-Virgo observing run. Accounting for selection biases,
we find the population of LIGO-Virgo black holes to be consistent with both homogeneity and
isotropy. We additionally find the data to constrain some directions of alignment more than others,
discuss the interpretation of this measurement and produce posteriors for the directions of total angular
momentum of all binaries in our set. While our current constraints are weak, the fact that such a small
number of detections can already yield a measurement suggests that this will be a powerful tool in
the future; we explore this prospect with a number of simulated catalogs of varying size.

[![arXiv](https://img.shields.io/badge/arXiv-2304.13254-crimson)](https://arxiv.org/abs/2304.13254)


## Installation
In order to compile the manuscript using `showyourwork`, you might need to install a custom backwards-compatible version of this package by doing:

```
pip install -U git+https://github.com/maxisi/showyourwork.git@gwisotropy
```

You can then simply build the article by calling the following from within the rootdir of this repository:

```
showyourwork build
```

For further details, see **[`src/scripts`](src/scripts)** and the _showyourwork_ documentation (linked below).

<p align="center">
<a href="https://github.com/showyourwork/showyourwork">
<img width = "450" src="https://raw.githubusercontent.com/showyourwork/.github/main/images/showyourwork.png" alt="showyourwork"/>
</a>
<br>
<br>
<a href="https://github.com/maxisi/gwisotropy/actions/workflows/build.yml">
<img src="https://github.com/maxisi/gwisotropy/actions/workflows/build.yml/badge.svg?branch=main" alt="Article status"/>
</a>
<a href="https://github.com/maxisi/gwisotropy/raw/main-pdf/arxiv.tar.gz">
<img src="https://img.shields.io/badge/article-tarball-blue.svg?style=flat" alt="Article tarball"/>
</a>
<a href="https://doi.org/10.5281/zenodo.7775266"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.7775266.svg" alt="DOI"></a>
<a href="https://github.com/maxisi/gwisotropy/raw/main-pdf/gwisotropy.pdf">
<img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
</a>
</p>

An open source scientific article created using the [showyourwork](https://github.com/showyourwork/showyourwork) workflow. Click [here](https://github.com/maxisi/gwisotropy/blob/main-pdf/gwisotropy.pdf) to view PDF without automatically downloading.
