import numpy as np
from astropy.cosmology import default_cosmology, z_at_value, Planck15
import lal
import lalsimulation as lalsim
from . import settings

RNG = np.random.default_rng()

# cosmo = default_cosmology.get_cosmology_from_string('Planck15')

# see https://docs.astropy.org/en/stable/api/astropy.cosmology.default_cosmology.html
with default_cosmology.set(Planck15):
    cosmo = default_cosmology.get()

# #############################################################################
# INTRINSIC
# #############################################################################

# global defaults
MMIN = settings.MMIN
MMAX = settings.MMAX
ZMIN = settings.ZMIN
ZMAX = settings.ZMAX

# See discussion below for justification of this 1/m1 choice
# https://wiki.ligo.org/CBC/RatesPop/O3aPopulationInjections#Will_F_proposal_38_discussion
def draw_m1q(m_min=MMIN, m_max=MMAX, rng=RNG):
    """ Draw m1 and q = m2/m1 from a distribution proportional to 1/m1 and
    uniform in q.
    """
    m1 = np.exp(rng.uniform(np.log(m_min), np.log(m_max)))
    q = rng.uniform()
    if m1*q < 5:
        m1, q = draw_m1q(m_min, m_max, rng)
    return m1, q

# draw spins from a sphere
def draw_iso_spins(n=1, rng=RNG):
    sx = rng.normal(size=n)
    sy = rng.normal(size=n)
    sz = rng.normal(size=n)
    s = np.sqrt(sx**2 + sy**2 + sz**2)
    new_s = rng.uniform(0, 0.9)
    return sx*new_s/s, sy*new_s/s, sz*new_s/s


# #############################################################################
# RATES
# #############################################################################

def md_sfr(z):
    """ Madau-Dickinson star formation rate. This rate is defined in the 
    source frame, so it has to be redshifted to translate into observed 
    quantities. It is defined as 1/(volume * time)
    """
    # TODO: shouln't this be 2.9 ?
    return (1+z)**2.7 / (1 + ((1+z)/(1+1.9))**5.6)

def z_grid(zmin=ZMIN, zmax=ZMAX):
    # define a grid of redshifts to integrate over
    zinterp = np.expm1(np.arange(np.log(zmin), np.log(1.0+zmax), 0.01))
    
    # set the cosmology
    
    # compute the volume element (per steradian) as a function of redshift for
    # each point in our redshift grid
    dVdz_interp = cosmo.differential_comoving_volume(zinterp)
    
    # compute number of sources per redshift bin and comoving time,
    # accounting for cosmological time dilation in the observed rate
    # (dV_c/dz) (d t_c / dt) = (dV_c/dz) / (1 + z)
    # i.e., events at higher redshifts are perceived to occur with increased 
    # delays, so that the observed rate is lower than the intrinsic one
    n = md_sfr(zinterp) * dVdz_interp / (1 + zinterp)
    
    # compute weights (z_weights) for redshift bins given the above
    # (we are more likely to detect sources coming from redshift bins
    # corresponding to more comoving volume and faster time)
    zwt = 0.5 * (n[:-1] + n[1:]) * np.diff(zinterp)
    zwt /= sum(zwt)
    # recenter redshift bins, from which we will draw given the above weights
    zdraw = 0.5 * (zinterp[:-1] + zinterp[1:])

    # for use further down, compute norm of z distribution
    znorm = np.trapz(4*np.pi*n, zinterp)
    return zdraw, zwt, znorm

def draw_z_sfr(zdraw=None, zwt=None, rng=RNG, **kws):
    """Draw redshifts according to the MDR distribution and the Vc weighting
    defined above for p(z).
    """
    if zdraw is None or zwt is None:
        zdraw, zwt, _ = z_grid(**kws)
    return rng.choice(zdraw, p=zwt)


# #############################################################################
# WEIGHTS
# #############################################################################

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

# #############################################################################
# SNR
# #############################################################################

def next_power_of_two(x):
    value = 1
    while value <= x:
        value = value  <<  1
    return value

def draw_snr(m1, q, dL_mpc, iota, psi, ra, dec, phi, s1x, s1y, s1z, s2x, s2y, s2z, 
             fstart=9., flow=10., fhigh=1024., fref=20., rng=RNG, asd_paths=None):
    """Function to draw waveform parameters based on fixed distributions
    and compute corresponding SNR at all detectors.
    
    Returns
    -------
    snrs: list
        list of ifo optimal SNRs
    """
    # compute secondary mass from m1 and q
    m2 = q*m1
    
    # draw random phase angle and hour angle
    # (we do this here rather than in the block below because we don't care 
    # about recording these parameters)
    gmst = 2*np.pi*rng.normal()
    
    # get expected waveform duration and corresponding delta_f
    T = next_power_of_two(lalsim.SimInspiralChirpTimeBound(fstart, m1*lal.MSUN_SI, m2*lal.MSUN_SI, 0.0, 0.0))
    dF = 1. / T

    # call waveform generator with long_asc_node = 0, which should be correct given how
    # XLALSimInspiralChooseFDWaveformFromCache calls XLALSimInspiralChooseFDWaveform in
    # https://docs.ligo.org/lscsoft/lalsuite/lalsimulation/group___l_a_l_sim_inspiral_waveform_cache__h.html#ga7cbe639459d24eb501f026a624e64ed3
    # These are the functions called by LALInference in 
    # https://git.ligo.org/lscsoft/lalsuite/-/blob/master/lalinference/lib/LALInferenceTemplate.c#L992
    # (presumably, bilby/RIFT do the same...)
    try:
        hp, hc = lalsim.SimInspiralChooseFDWaveform(m1*lal.MSUN_SI, m2*lal.MSUN_SI,
                                                    s1x, s1y, s1z,
                                                    s2x, s2y, s2z,
                                                    float(dL_mpc*1e6*lal.PC_SI), iota,
                                                    phi, 0.0, 0.0, 0.0, dF,
                                                    fstart, fhigh, fref,
                                                    lal.CreateDict(), lalsim.IMRPhenomPv3)
    except RuntimeError:
        print(f"Parameters failed: {m1}, {m2}, {s1x}, {s1y}, {s1z}, {s2x}, {s2y}, {s2z}, {dL_mpc}, {iota}, {phi}")
        return [0]*len(asd_paths)

    # compute IFO optimal SNRs for the drawn signal
    SNRs = []
    for ifo, asd_path in asd_paths.items():
        d = lal.cached_detector_by_prefix[ifo]
        Fp, Fc = lal.ComputeDetAMResponse(d.response, ra, dec, psi, gmst)
        h = lal.CreateCOMPLEX16FrequencySeries("h in detector", lal.LIGOTimeGPS(), 
                                               fstart, dF, hp.sampleUnits, hp.data.length)
        h.data.data = Fp * hp.data.data + Fc*hc.data.data

        # create PSD from ASD file 
        # https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/group___l_a_l_sim_noise_p_s_d__c.html#ga67d250556e4e8647c50d379364eb7911
        psd = lal.CreateREAL8FrequencySeries("detector PSD", lal.LIGOTimeGPS(),
                                             fstart, dF, lal.DimensionlessUnit,
                                             hp.data.length)
        lalsim.SimNoisePSDFromFile(psd, fstart, asd_path)
        SNRs.append(lalsim.MeasureSNRFD(h, psd, flow, fhigh))
    return SNRs
