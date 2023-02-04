import numpy as np
from .settings import SIGMA_RAW

def draw_prior(sigma=SIGMA_RAW, ndraw=100000, rng=None, seed=None):
    """Draw `ndraw` Cartesian vectors from the prior.

    Arguments
    ---------
    sigma : float
        prior standard deviation, defaults to 0.4
    ndraw : int
        number of vectors to draw, def. 100000
    rng : np.random.Generator
        random number generator, def. `np.random.default_rng`
    seed : int
        seed for `np.random.default_rng` if `rng` is not provided, def. None

    Returns
    -------
    v : np.array
        array of shape `(ndraw, 3)` containing the Cartesian vectors drawn from
        the prior.
    """

    if rng is None:
        rng = np.random.default_rngs(seed)

    # we start with a Gaussian prior on the components of vL_raw
    v_raw = np.random.normal(0, sigma, size=(ndraw, 3))

    # then, we can compute be and its norm from the definition above
    v = (v_raw.T / np.sqrt(1 + np.einsum("ij,ij->i", v_raw, v_raw))).T
    return v

