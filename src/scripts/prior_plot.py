import paths
import seaborn as sns
from utils.prior import draw_prior
import numpy as np

sns.set(context='notebook', palette='colorblind')

RNG = np.random.default_rng(12345)

# draw vectors from prior
v = draw_prior(ndraw=100000, rng=RNG)

# plot 2D slice
with sns.axes_style("ticks"):
    g = sns.jointplot(x=v[:,0], y=v[:,1], kind='hex');
    g.ax_joint.set_xlabel(r"$\vec{v}_x$")
    g.ax_joint.set_ylabel(r"$\vec{v}_y$");

    p = paths.figures / "prior.pdf"
    g.fig.savefig(p, bbox_inches='tight')
    print(f"Saved: {p}")
