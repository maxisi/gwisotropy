# -*- coding: utf-8 -*-
#
#       Copyright 2023
#       Maximiliano Isi <max.isi@ligo.org>
#       Will M. Farr <will.farr@ligo.org>
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

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

