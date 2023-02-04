import paths
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import arviz as az
import utils
import utils.plots
from utils.prior import draw_prior

sns.set(context='notebook', palette='colorblind')

RNG = np.random.default_rng(12345)


###############################################################################
# LOAD FIT
###############################################################################

fit = az.from_netcdf(paths.result)

x = fit.posterior.vN.values
nhats = x.reshape(np.prod(x.shape[:2]), 3)

x = fit.posterior.vL.values
jhats = x.reshape(np.prod(x.shape[:2]), 3)

# ----------------------------------------------------------------------------
# draw vectors from prior

vprior = draw_prior(ndraw=len(jhats), rng=RNG)

###############################################################################
# PLOT
###############################################################################

lkws = dict(c='gray', ls='--')

with sns.axes_style("ticks"):
    fig = plt.figure(figsize=utils.plots.figsize_column)
    sns.kdeplot(np.linalg.norm(jhats, axis=1), label="$J$ posterior", lw=3)
    sns.kdeplot(np.linalg.norm(nhats, axis=1), label="$N$ posterior", lw=3)
    sns.kdeplot(np.linalg.norm(vprior, axis=1), label="prior", color='0.8',
                lw=3, ls='--', zorder=-100);
    plt.xlabel(r"$|\vec{v}_{J/N}|$");
    plt.legend();

    p = paths.figures / f"jn_norm.pdf"
    fig.savefig(p, bbox_inches='tight')
    print(f"Saved: {p}")
