import paths
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import arviz as az
from scipy.stats import gaussian_kde
import utils
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

vdict = {'n': nhats, 'j': jhats, 'prior': vprior}

###############################################################################
# PLOT
###############################################################################

lkws = dict(c='gray', ls='--')

for k, vs in vdict.items():
    # KDE vectors to later color points by density
    x ,y ,z  = vs.T
    kde = gaussian_kde([x, y, z])
    vcs = kde([x, y, z])
    
    with sns.axes_style("whitegrid", {"grid.linestyle": ':'}):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
    
        ax.plot3D([0, 0], [0,0], [-1,1], **lkws)
        ax.plot3D([0, 0], [-1,1], [0,0], **lkws)
        ax.plot3D([-1,1], [0, 0], [0,0], **lkws)
    
        ax.scatter(x, y, z, c=vcs, cmap='viridis', alpha=0.3, marker='.', lw=0,
                   rasterized=True)
        ax.scatter([0], [0], [0], c='k', marker='P', s=40)
    
        for i in 'xyz':
            getattr(ax, f'set_{i}lim')(-1,1)
            getattr(ax, f'set_{i}ticks')([-1,0,1])
            getattr(ax, f'set_{i}label')(r"$%s$" % i, labelpad=-1)
        ax.set_zlabel("$z$", labelpad=-8)
        ax.tick_params(axis='both', which='major', pad=-1)

    p = paths.figures / f"density_3d_{k}.pdf"
    fig.savefig(p, bbox_inches='tight', dpi=600)
    print(f"Saved: {p}")
