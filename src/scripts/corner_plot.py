import paths
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import utils
import utils.plots
import pandas as pd
import arviz as az
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.lines import Line2D

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

###############################################################################
# PLOT
###############################################################################

ckeys = [r'$v_{N,x}$', r'$v_{N,y}$', r'$v_{N,z}$',
         r'$v_{J,x}$', r'$v_{J,y}$', r'$v_{J,z}$']
df = pd.DataFrame(np.hstack((nhats, jhats)), columns=ckeys)

lkws = dict(ls=':', c='gray', alpha=0.5)
with sns.axes_style("ticks"):
    pg = sns.PairGrid(df, corner=True)
    pg.map_diag(sns.kdeplot)
    pg.map_lower(sns.kdeplot, fill=False)
    for i, axs in enumerate(pg.axes):
        for j, ax in enumerate(axs):
            if ax:
                ax.set_xlim(-1, 1)
                ax.axvline(0, **lkws)
                if i != j:
                    ax.axhline(0, **lkws)

sky_kws = dict(cmap='viridis', s=20, rasterized=True)
akws = dict(arrowstyle="->", color=sns.color_palette()[0], alpha=0.3, lw=2.5)
with sns.axes_style("whitegrid", {"grid.linestyle": ':'}):
    axins = pg.axes[1,1].inset_axes([3.25, 0.5, 2.5, 1.5], projection='mollweide')
    utils.plots.sky_scatter(fit.posterior.vN.values, ax=axins, **sky_kws)
    utils.plots.add_colorbar(axins, key='N', cmap=sky_kws['cmap'])
    axins.annotate("", xy=(-0.75, 0.5), xytext=(-0.25, 0.5),
                   arrowprops=akws, xycoords='axes fraction')

    axins = pg.axes[4,4].inset_axes([-0.25, 2, 2.5, 1.5], projection='mollweide')
    utils.plots.sky_scatter(fit.posterior.vL.values, ax=axins, **sky_kws)
    utils.plots.add_colorbar(axins, key='J', cmap=sky_kws['cmap']);
    axins.annotate("", xy=(0.5, -0.75), xytext=(0.5, -0.25),
               arrowprops=akws, xycoords='axes fraction')

xys = [([-0.02, -0.02], [0.52, 0.95]), ([0.52, 1], [0, 0])]
lines = [Line2D(x, y, lw=5., alpha=0.3, figure=pg.fig, transform=pg.fig.transFigure,
               ) for x, y in xys]
pg.fig.patches.extend(lines)

p = paths.figures / "jn_corner.pdf"
pg.fig.savefig(p, bbox_inches="tight")
print(f"Saved: {p}")