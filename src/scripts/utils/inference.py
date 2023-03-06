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

import pymc as pm
import numpy as np
import aesara.tensor as at
import arviz as az
import scipy.stats as ss
from . import settings

def make_model(location_posterior_stack, orientation_posterior_stack,
               location_selection_stack, orientation_selection_stack,
               pdraw, Ndraw, sigma_raw=settings.SIGMA_RAW):
    """Construct isotropy probability model in PyMC.

    The likelihood function takes the form:

    . math::

        \ln \mathcal{L} \sim  \frac{1}{16 \pi^2} \left(1 +
        \vec{v}_n\cdot\vec{n}\right) \left(1 + \vec{v}_{L}\cdot\vec{L}\right)

    where :math:`\vec{v}_n` and :math:`\vec{v}_L` are 3D vectors defining two
    special directions in galactic coordinates, respectively defining a
    preferred sky location (:math:`\hat{n}`) or orientation (:math:`\hat{L}`);
    here the orientation is defined in terms of the _orbital_ angular momentum
    :math:`\vec{L}`, but could equivalently be replaced by the _total_ angular
    momentum :math:`\vec{J}`.
 
    Isotropy is recovered for :math:`\vec{v}_n = \vec{v}_{L,J} = 0`.

    This is a model for the intrinsic distribution sources, not detected
    distribution. To realize this, we take into account the selection function
    encoding the detectability of different parameters as a result of imperfect
    detector sensitivity.

    In order to account for this, we must add a term to the log likelihood like

    . math::
        \ln \mathcal{L}_\mathrm{sel} = - N_\mathrm{obs} \ln \alpha(\vec{v}_n,
        \vec{v}_{L,J})

    where :math:`N_\mathrm{obs}` is the number of observed events and
    :math:`\alpha(\vec{v}_n, \vec{v}_{L,J})` is the selection term, estimated
    from the injection set by

    . math::
        \hat{\alpha} = \frac{1}{N_{\rm draw}} \sum_{i=1}^{N_{\rm det}}
        \frac{p(\hat{n}, \hat{L} \mid \vec{v}_n,
        \vec{v}_{L,J})}{p_\mathrm{draw}(\hat{n}, \hat{L})}

    where the sum is over the :math:`N_\mathrm{det}` simulated signals that
    were detected out of the original :math:`N_\mathrm{draw}` injections in the
    set, and :math:`p_\mathrm{draw}` is the original probability of drawing a
    given set of parameters from the injection set in the first place.

    See https://arxiv.org/abs/2204.00461 , https://arxiv.org/abs/1904.10879

    Arguments
    ---------
    location_posterior_stack : array
        three-dimensional array with Cartesian vectors for the sky location
        posteriors of each event, shaped as ``(n_events, n_samples, 3)``.
    orientation_posterior_stack : array
        same as ``location_posterior_stack`` but for the vectors encoding the 
        binary orientation (total or orbital angular momentum).
    location_selection_stack : array
        two-dimensional array with Cartesian vectors for the sky location
        of each detected event from the injection campaign to quantifyg 
        selection effects; shaped as ``(n_det, 3)``.
    orientation_selection_stack : array
        same as ``location_selection_stack`` but for the vectors encoding the 
        binary orientation (total or orbital angular momentum).
    pdraw : array
        draw probability for each detected signal in the selection set.
    ndraw : int
        number of total injections drawn (i.e., the sum of detected and
        nondetected signals).
    sigma_raw : float
        standard deviation of Gaussian prior for "raw" dipole vectors ``v_raw``
        used in sampling.
    """
    with pm.Model() as model:
        Nobs = len(location_posterior_stack)

        # draw random special location vector, rescale
        vN_raw = pm.Normal('vN_raw', mu=0, sigma=sigma_raw, shape=3)
        vN = vN_raw / at.sqrt(1 + at.dot(vN_raw, vN_raw))

        # draw random special orientation vector, rescale
        vL_raw = pm.Normal('vL_raw', mu=0, sigma=sigma_raw, shape=3)
        vL = vL_raw / at.sqrt(1 + at.dot(vL_raw, vL_raw))

        # compute spherical density for both location and orientation
        def sph_density(Ns, Ls, axes=(0,2)):
            n_ip = at.tensordot(vN, Ns, axes)
            l_ip = at.tensordot(vL, Ls, axes)
            return 1/(16*np.pi**2)*(1 + n_ip)*(1 + l_ip)

        # compute likelihood given observed posteriors
        evt_wts = sph_density(location_posterior_stack,
                              orientation_posterior_stack)
        evt_log_mean_wts = at.log(at.mean(evt_wts, 1))

        # manually add totally log-likelihood to the potential
        pm.Potential('evt_wts_lnlike', at.sum(evt_log_mean_wts))

        # evaluate the selection function (mu_sel is `alpha`), and estimate
        # the accuracy of the approximation (through var_sel)
        # see https://arxiv.org/abs/1904.10879
        sel_wts = sph_density(location_selection_stack,
                              orientation_selection_stack,
                              axes=(0,1))
        mu_sel = at.sum(sel_wts / pdraw) / Ndraw
        var_sel = at.sum(sel_wts**2 / pdraw**2) / Ndraw**2 - mu_sel**2 / Ndraw
        Neff_sel = mu_sel**2 / var_sel

        # manually add selection weight to the potential
        pm.Potential('sel_wts_lnlike', -Nobs*at.log(mu_sel))

        # generated quantities
        pm.Deterministic("vL", vL)
        pm.Deterministic("vN", vN)
        pm.Deterministic("Neff_sel", Neff_sel)
    return model

def draw_prior(sigma=1/np.sqrt(5), ndraw=100000, rng=None, seed=None):
    """Draw `ndraw` Cartesian vectors from the prior.

    Arguments
    ---------
    sigma : float
        prior standard deviation, defaults to 1/sqrt(5).
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

def cl_origin(vecs):
    """Compute the credibile level of the origin, defined as the fraction of
    samples with probability density greater than the origin.
    
    Arguments
    ---------
    vecs : array
        Nx3 array containing three Cartesian vector components for N samples.

    Returns
    -------
    cl : float
        credible level.
    """
    kde = ss.gaussian_kde(vecs.T)
    po = kde([0,0,0])
    ps = kde(vecs.T)
    return np.count_nonzero(ps > po) / len(ps)
