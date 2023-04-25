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
from astropy.cosmology import default_cosmology, z_at_value, Planck15
from . import settings
import scipy.stats as ss

RNG = np.random.default_rng()

# see https://docs.astropy.org/en/stable/api/astropy.cosmology.default_cosmology.html
with default_cosmology.set(Planck15):
    cosmo = default_cosmology.get()

# global defaults
MMIN = settings.MMIN
MMAX = settings.MMAX
ZMIN = settings.ZMIN
ZMAX = settings.ZMAX

# #############################################################################
# WEIGHTS
# #############################################################################

def md_sfr(z):
    """ Madau-Dickinson star formation rate. This rate is defined in the 
    source frame, so it has to be redshifted to translate into observed 
    quantities. It is defined as 1/(volume * time)
    """
    return (1+z)**2.7 / (1 + ((1+z)/(1+1.9))**5.6)

def li_prior_wt(m1_src, q, z, ra, sin_dec, cos_iota, psi, cosmo_prior=False):
    """ Returns the prior weight applied by LALInference in the given coordinate
    system.
    
    If `from_cosmo = True`, assumes we start from samples that have been 
    reweighted in redshift to be uniform in comoving volume, so it only
    applies the mass reweighting.
    
    The LALInference prior function in the source-frame masses and redshift
    is just given by the following Jacobian (see Appendix C in [Abbott et al]
    (https://iopscience.iop.org/article/10.3847/2041-8213/ab3800/pdf)):

    . math::
    p(m_1^{src}, m_2^{src}, z) = p(m_1, m_2, d_L) 
    \left|\frac{\partial(m_1, m_2, d_L)}
    {\partial (m_1^{src}, m_2^{src}, z)}\right| \propto d^2_L(z) 
    \frac{\partial m_1}{\partial m_1^{src}} 
    \frac{\partial m_2}{\partial m_2^{src}} 
    \frac{\partial d_L}{\partial z}
    = d^2_L(z) (1+z)^2 \frac{\partial d_L}{\partial z}

    because :math:`m_{1/2} = (1+z) m_{1/2}^{src}`. Now, to compute the
    :math:`\partial d_L/ \partial z` term, we need knowledge of cosmology.
    [Hogg (1999)](https://arxiv.org/abs/astro-ph/9905116) shows that,

    . math::
    \frac{\partial d_L}{\partial z} = \frac{d_L}{1+z} + \left(1 + z\right) 
    \frac{d_H}{E(z)}

    where :math:`d_H = c/H_0 and 
    :math:`E(z) = \sqrt{\Omega_M (1 + z)^3 + \Omega_\Lambda}`.

    We actually use :math:`q` and :math:`m_1`, instead of :math:`m_1` and
    :math:`m_2`, so the Jacobian we apply is:

    . math::
    p(m_1^{src}, q, z) = p(m_1, q, d_L) \left|\frac{\partial(m_1, q, d_L)}
    {\partial (m_1^{src}, q, z)}\right| = p(m_1, m_2, d_L)
    \left|\frac{\partial m_2}{\partial q}\right| 
    \left|\frac{\partial(m_1, q, d_L)}{\partial (m_1^{src}, q, z)}\right|
    \propto d^2_L(z) m_1 \frac{\partial m_1}{\partial m_1^{src}} 
    \frac{\partial d_L}{\partial z}
    = d^2_L(z)\, m_1\, (1+z) \frac{\partial d_L}{\partial z} =
    d^2_L(z)\, m_1^{src}\, (1+z)^2 \frac{\partial d_L}{\partial z} 
    """
    if cosmo_prior:
        dVdz = 4*np.pi*cosmo.differential_comoving_volume(z).value
        z_wt = dVdz/(1+z)
    else:
        dL = cosmo.luminosity_distance(z)
        ddL_dz = dL/(1+z) + (1+z)*cosmo.hubble_distance / cosmo.efunc(z)
        z_wt = dL**2 * ddL_dz
    mass_wt = (1+z)**2*m1_src
    return mass_wt*z_wt

def mass_redshift_pop_wt(m1_src, q, z, m_min=MMIN, m_max=MMAX, z_max=ZMAX):
    """Returns the weighting for the default population (also the injected population)
    in mass and redshift: :math:`m_{min} \\, M_\\odot < m_2 < m_1 < m_{max} \\, M_\\odot`,
    flat in log-m1, flat in q, and merger rate proportional to the [Madau & Dickinson
    (2014)](https://ui.adsabs.harvard.edu/abs/2014ARA%26A..52..415M/abstract) SFR.
    
    If `from_cosmo = True`, assumes we start from samples that have been 
    reweighted in redshift to be uniform in comoving volume, so it only
    applies the mass reweighting.
    """
    m1_src = np.atleast_1d(m1_src)
    q = np.atleast_1d(q)
    z = np.atleast_1d(z)
    mask = (m1_src > m_max) | (m1_src*q < m_min) | (z > z_max) | (z < 0) | (q < 0) | (q > 1)
    md_wt = md_sfr(z)/m1_src
    dVdz = np.zeros_like(z)
    dVdz[~mask] = 4*np.pi*cosmo.differential_comoving_volume(z[~mask]).value
    z_wt = dVdz/(1+z)
    return md_wt * z_wt

# EXTRA: R&P DISTRIBUTIONS

def rp_chi_wt(chi1, chi2, refdf):
    # a Beta distribution, Eq. (B19) in arXiv:2111.03634
    spin_dist = ss.beta(refdf['alpha_chi'], refdf['beta_chi'])
    return spin_dist.pdf(chi1)*spin_dist.pdf(chi2)

def rp_tilt_wt(costh1, costh2, refdf):
    # Eq. (B20) in arXiv:2111.03634
    xi = refdf['xi_spin']
    sigma = refdf['sigma_spin']
    myclip_a, myclip_b = -1, 1
    a, b = myclip_a/sigma, myclip_b/sigma
    truncnorm = ss.truncnorm(a=a, b=b, scale=sigma)
    isotropic = 1/4.
    return xi*truncnorm.pdf(costh1)*truncnorm.pdf(costh2) + (1-xi)*isotropic

def rp_wt(m1_src, q, z, chi1, chi2, costh1, costh2, refdf=None, m_min=MMIN,
          m_max=MMAX, z_max=ZMAX):
    spin_wt = rp_chi_wt(chi1, chi2, refdf)*rp_tilt_wt(costh1, costh2, refdf)
    
    m1_src = np.atleast_1d(m1_src)
    q = np.atleast_1d(q)
    z = np.atleast_1d(z)
    mask = (m1_src > m_max) | (m1_src*q < m_min) | (z > z_max) | (z < 0) |\
           (q < 0) | (q > 1)
    md_wt = md_sfr(z)/m1_src
    dVdz = np.zeros_like(z)
    dVdz[~mask] = 4*np.pi*cosmo.differential_comoving_volume(z[~mask]).value
    z_wt = dVdz/(1+z)

    return spin_wt * md_wt * z_wt

# def smoothing_f(mp, delta_m):
# 		# Eq. (B6) in arXiv:2111.03634
#     return np.exp(delta_m/mp + delta_m/(mp-delta_m))
# 
# def smoothing_s(m, mmin, delta_m):
# 		# Eq. (B5) in arXiv:2111.03634
#     s = np.ones_like(m)
#     s[m < mmin] = 0.
#     mask = (mmin <= m) & (m < mmin + delta_m)
#     s[mask] = 1./(smoothing_f(m[mask] - mmin, delta_m) + 1)
#     return s 
# 
# def rp_mass_wt(m1_src, q, refdf):
#     # Eq. (B4) in arXiv:2111.03634
#     # MASS
#     m1_src = np.atleast_1d(m1_src)
#     q = np.atleast_1d(q)
#     alpha, mmax, mmin, delta_m = refdf['alpha'], refdf['mmax'],\
#                                  refdf['mmin'], refdf['delta_m']
#     # power law
#     power_law_norm = (mmin**(1-alpha) - mmax**(1-alpha))/(alpha-1)
#     power_law = m1_src**(-alpha) / power_law_norm
#     power_law *= np.heaviside(m1_src - mmax, 1)*(1 - refdf['lam'])
#     # Gaussian 
#     gaussian = ss.norm.pdf(m1_src, refdf['mpp'], refdf['sigpp'])
#     gaussian *= refdf['lam']
#     p_m1_src = (power_law + gaussian)*smoothing_s(m1_src, mmin, delta_m)
# 
#     # MASS RATIO
#     p_q = q**refdf['beta'] * smoothing_s(q*m1_src, mmin, delta_m) * (1 + refdf['beta'])
# 
#     return p_m1_src * p_q

