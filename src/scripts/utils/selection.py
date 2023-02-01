import numpy as np
import pandas as pd
from . import pops
import os
from tqdm import tqdm
from astropy.cosmology import default_cosmology, z_at_value, Planck15
# set the cosmology
with default_cosmology.set(Planck15):
    cosmo = default_cosmology.get()

def save_selection_data(df, Ndraw, path):
    df.to_hdf(path, 'table')
    with pd.HDFStore(path) as store:
        store.get_storer('table').attrs.metadata = {'Ndraw': Ndraw}

def load_selection_data(path):
    df = pd.read_hdf(path, 'table')
    with pd.HDFStore(path) as store:
        Ndraw = store.get_storer('table').attrs.metadata['Ndraw']
    return df, Ndraw

RNG = np.random.default_rng()
def draw_signals(Ndetected=1000, snr_thresh=9.0, path=None, m_min=3., m_max=150.,
                 zdraw=None, zwt=None, z_min=None, z_max=None, rng=RNG, **kws):
    # # requested number of detected sources
    # Ndetected = 1000
    # # SNR threshold defining detection
    # snr_thresh = 9.0

    # these are the parameters we will record for each detected injection
    row_keys = ['m1', 'q', 'z', 'cosiota', 'psi', 'ra', 'sindec', 'phase',
                'spin_1x', 'spin_1y', 'spin_1z', 'spin_2x', 'spin_2y', 'spin_2z',
                'pdraw', 'pdrawangle']
    
    if path is not None and os.path.isfile(path):
        df_sel, Ndraw = load_selection_data(path)
    else:
        Ndraw = 0
        with tqdm(total=Ndetected) as pbar:
            if zdraw is None or zwt is None:
                zdraw, zwt, znorm = pops.z_grid(z_min, z_max)

            # iterate until we reach the desired number of detections
            rows = []
            while len(rows) < Ndetected:
                Ndraw += 1
    
                # draw intrinsic parameters
                m1, q = pops.draw_m1q(m_min, m_max, rng)
                
                # draw redshfit and compute luminosity distance
                z = pops.draw_z_sfr(zdraw, zwt, rng=rng)
                d = cosmo.luminosity_distance(z)
                
                # draw spins
                s1x, s1y, s1z = pops.draw_iso_spins()
                s2x, s2y, s2z = pops.draw_iso_spins()
                s1x = s1x[0]
                s1y = s1y[0]
                s1z = s1z[0]
                s2x = s2x[0]
                s2y = s2y[0]
                s2z = s2z[0]
    
                # draw orientation parameters and phase
                cos_iota = rng.uniform(-1, 1)
                psi = rng.uniform(0, 2*np.pi)
                ra = rng.uniform(0, 2*np.pi)
                sin_dec = rng.uniform(-1, 1)
                phi = rng.uniform(0, 2*np.pi)
    
                iota = np.arccos(cos_iota)
                dec = np.arcsin(sin_dec)
    
                # compute (optimal) individual detector SNRs
                # (will generate a random phase and hour angle)
                snrs = np.array(pops.draw_snr(m1*(1+z), q, d.value, iota, psi, ra, dec,
                                              phi, s1x, s1y, s1z, s2x, s2y, s2z, rng=rng, **kws))
                
                # compute network SNR, adding some noise to simulate MF SNR
                # instead of optimal
                snr = np.linalg.norm(snrs + rng.normal(size=len(snrs)))
    
                if snr > snr_thresh:
                    # the signal was detected, so let's compute and record
                    # its "prior" (according to the injection distribution; aka, 
                    # the draw probability)
    
                    # the mass probability densities are given by draw_m1q()
                    # replicate that here: uniform in log(m1) and uniform in q
                    pm = 1/(m1*(np.log(m_max)-np.log(m_min)))
                    pq = 1
                    
                    # the redshift PDF follows MD in the comoving space, 
                    # according to draw_z_sfr()
                    dVdz = cosmo.differential_comoving_volume(z)
                    pz = pops.md_sfr(z)*dVdz/(1 + z)/znorm
                    
                    # the distributions on the angles is isotropic
                    pcosiota = 1/2
                    ppsi = 1/(2*np.pi)
                    pra = 1/(2*np.pi)
                    psindec = 1/2
                    pdraw_angle = pcosiota*ppsi*pra*psindec
                    
                    # the distribution on the phase is also isotropic
                    pphi = 1/(2*np.pi)
                    
                    # the total probability is the product of all
                    pdraw = pm*pq*pz*pdraw_angle*pphi
    
                    row = dict(zip(row_keys, [m1, q, z, cos_iota, psi, ra, sin_dec, phi,
                                              s1x, s1y, s1z, s2x, s2y, s2z,
                                              pdraw, pdraw_angle]))
                    rows.append(row)
    
                    pbar.update(1)
        df_sel = pd.DataFrame(rows)
        df_sel['pdraw'] = df_sel['pdraw'].map(lambda x:  x.value)
        if path is not None:
            save_selection_data(df_sel, Ndraw, path)
    return df_sel, Ndraw
