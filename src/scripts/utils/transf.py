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
from tqdm import tqdm
import pandas as pd

MTSUN_SI = 4.925491025543576248120047206269e-06

def fix_psi(samples, duplicate=False, rng=None):
    """Restore polarization angle to full range [0, 2*pi].

    Arguments
    ---------
    samples: DataFrame, ndarray
        original sampled
    duplicate: bool
        whether to fix by doubling the number of samples (def., False).
    rng:
        random number generator (opt.)

    Returns
    -------
    new_samples: DataFrame, ndarray
        copy of `samples` with fixed psi.
    """
    if rng is None:
        rng = np.random.default_rng()
    new_samples = samples.copy()
    new_samples['psi'] = new_samples['psi'] % np.pi
    if duplicate:
        # create copy of samples with psi -> psi + pi
        samples_copy = new_samples.copy()
        samples_copy['psi'] = samples_copy['psi'] + np.pi
        new_samples = np.hstack([new_samples, samples_copy])
        rng.shuffle(new_samples)
    else:
        # randomly add pi or 0 to each sample
        new_samples['psi'] += rng.integers(0, 2, len(new_samples))*np.pi
    return new_samples


def rotation_matrix(angle, direction):
    """Return the matrix effecting a rotation by angle, counterclockwise
    around the direction n (i.e. in the sense of the right-hand rule).
    """
    n = direction / np.linalg.norm(direction)
    costh, sinth = np.cos(angle), np.sin(angle)
    p = 1 - costh
    R = np.array([[costh + n[0]*n[0]*p,      n[0]*n[1]*p - n[2]*sinth, n[0]*n[2]*p + n[1]*sinth],
                [n[0]*n[1]*p + n[2]*sinth, costh + n[1]*n[1]*p,      n[1]*n[2]*p - n[0]*sinth],
                [n[0]*n[2]*p - n[1]*sinth, n[2]*n[1]*p + n[0]*sinth, costh + n[2]*n[2]*p]])
    return R
    
def get_location_vector(alpha, sindelta, *args):
    """Celestial coordinates are defined with alpha and delta as azimuthal 
    and (co)polar angles respectively, such that alpha = 0 corresponds to a
    source aligned with the March (vernal) equinox, which serves as X axis;
    delta = 0 corresponds to a source at the celestial equator, and so 
    sin(delta) is the Z coordinate along the celestial North
    """
    cosdelta = np.sqrt(1 - sindelta**2)
    x = cosdelta*np.cos(alpha)
    y = cosdelta*np.sin(alpha)
    z = sindelta
    return np.array([x, y, z])
    
def get_orientation_vector(alpha, sindelta, cosiota, psi, full_output=False):
    """iota is the polar angle between the line of sight and the orbital
    angular momentum, whereas psi is the in-sky angle between the angular
    momentum and the celestial East
    """
    # define celestial north
    north = np.array([0, 0, 1])
    # get wave-vector, which is Z in the waveframe
    n = get_location_vector(alpha, sindelta)
    k = -n
    # get local celestial West
    if abs(sindelta) == 1:
        # source at celestial North, need to disambiguate
        # waveframe orientation
        west = np.array([-np.cos(alpha+np.pi/2), -np.sin(alpha+np.pi/2), 0])
    else:
        west = np.cross(north, k)
        west /= np.linalg.norm(west)
    # the waveframe X is this vector rotated counterclockwise by psi around k,
    # as seen from Earth (i.e., following the right hand rule with thumb 
    # pointing towards k)
    wx = np.dot(rotation_matrix(psi, k), west)
    # finally rotate this vector, which lives in the plane of the sky,
    # such that it lies at an angle iota from z, i.e. rotate around Y
    wy = np.cross(k, wx)
    L = np.dot(rotation_matrix(np.arccos(cosiota), wy), k)
    if full_output:
        return L, west, wx/np.linalg.norm(wx), wy/np.linalg.norm(wy)
    else:
        return L

def SimInspiralLN(M, eta, v):
    # https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/_l_a_l_sim_inspiral_p_n_coefficients_8c_source.html#l02173
    return M*M*eta/v

def SimInspiralL_2PN(eta):
    # https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/_l_a_l_sim_inspiral_p_n_coefficients_8c_source.html#l02181
    return 1.5 + eta/6.;

def get_jhat_components(m1, m2, spin1_x, spin1_y, spin1_z, spin2_x, spin2_y, spin2_z, f_ref):
    """Get the components of J in a frame with L as the z axis, and with y
    from body 2 to body 1. This is the frame in which component spins are
    specified in the LALSimulation convention---see left panel of Fig. 1 
    in https://arxiv.org/pdf/2106.06492.pdf

    See the following LALSimulation function for clarification
    https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/_l_a_l_sim_inspiral_8c_source.html#l05885
    
    Arguments
    ---------
    m1 : float
        first component mass in solar masses
    m2 : float
        second component mass in solar masses
    spin{1,2}_{x,y,z} : float
        dimensionless spin components
    f_ref: float
        reference frequency in Hz
    """
    # the spin components are given in a frame where L points along Z and
    # X along the line from body 2 to body 1. Let's get J in that frame
    # To do that, we need the magnitude of L
    
    # get v parameter at reference point
    # see https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/_l_a_l_sim_inspiral_8c_source.html#l05923
    m = m1 + m2
    v0 = np.cbrt( m * MTSUN_SI * np.pi * f_ref );
    
    eta = m1*m2/m/m;

    Lmag = SimInspiralLN(m, eta, v0)*(1. + v0*v0*SimInspiralL_2PN(eta));
    
    # We also need the magnitudes of the dimensionfull spins
    s1x = m1 * m1 * spin1_x;
    s1y = m1 * m1 * spin1_y;
    s1z = m1 * m1 * spin1_z;
    s2x = m2 * m2 * spin2_x;
    s2y = m2 * m2 * spin2_y;
    s2z = m2 * m2 * spin2_z;
    
    # Now we are almost done! J is trivial to compute
    Jx = s1x + s2x;
    Jy = s1y + s2y;
    Jz = Lmag + s1z + s2z;
    
    # Let's normalize it
    Jnorm = np.sqrt( Jx*Jx + Jy*Jy + Jz*Jz);
    Jhatx = Jx / Jnorm;
    Jhaty = Jy / Jnorm;
    Jhatz = Jz / Jnorm;

    return Jhatx, Jhaty, Jhatz

def get_total_orientation_vector(alpha, sindelta, cosiota, psi, m1, m2, 
                                 spin1_x, spin1_y, spin1_z,
                                 spin2_x, spin2_y, spin2_z,
                                 phi_ref, f_ref, full_output=False):
    """Return the direction of the total angular momentum
    
    iota is the polar angle between the line of sight and the orbital
    angular momentum, whereas psi is the in-sky angle between the angular
    momentum and the celestial East
    """
    # get direction of orbital angular momentum
    L, _, _, wy = get_orientation_vector(alpha, sindelta, cosiota, psi,
                                         full_output=True)
    # get components of J in L-frame
    Jhatx, Jhaty, Jhatz = get_jhat_components(m1, m2, spin1_x, spin1_y, spin1_z,
                                              spin2_x, spin2_y, spin2_z, f_ref)
    # the Jhatz component points along L;
    # the Jhatx component points along a line from body 2 to body 1 in the 
    # orbital plane, which subtends  an angle phiRef from the line of nodes, 
    # i.e., x axis in
    # https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/group__lalsimulation__inspiral.html
    # Therefore, we get xhat by rotating the line of nodes (which is wy in the
    # LALSuite convention because Omega = pi/2) by phi_ref around L
    xhat = np.dot(rotation_matrix(phi_ref, L), wy)
    xhat /= np.linalg.norm(xhat)
    
    # the Jhaty component points along yhat, which is the vector that completes
    # the (xhat, yhat, L) triad
    yhat = np.cross(L, xhat)
    yhat /= np.linalg.norm(yhat)
    
    J = Jhatx[...,np.newaxis]*xhat[np.newaxis,...] + \
        Jhaty[...,np.newaxis]*yhat[np.newaxis,...] + \
        Jhatz[...,np.newaxis]*L[np.newaxis,...]/np.linalg.norm(L)

    return np.squeeze(J)


ks1 = ['ra', 'dec', 'iota', 'psi', 'mass_1', 'mass_2', 
       'spin_1x', 'spin_1y', 'spin_1z',
       'spin_2x', 'spin_2y', 'spin_2z',
       'phase',]
ks2 = ['ra', 'sindec', 'cosiota', 'psi']
ks3 = ['ra', 'sindec', 'cosiota', 'psi', 'mass_1', 'mass_2', 
       'spin_1x', 'spin_1y', 'spin_1z',
       'spin_2x', 'spin_2y', 'spin_2z',
       'phase', 'fref']

def get_vectors(rwt_samples_dict, fref=20.):
    """Get location (N) and orientation (J, L) vectors in Cartesian celestial 
    coordinates for a number of samples drawn in parameter estimation for each
    event.
    
    Arguments
    ---------
    rwt_samples_dict : dict
        dictionary of samples per event, like `{event: samples}`, where
        `samples` is a `pd.DataFrame` containing PE samples for `event`.
    fref : float
        reference frequency for spin orientations.
        
    Returns
    -------
    vecs_location : list
        binary locations; list of tuples with entries of dimension `(Nsamp, 3)`, 
        where `Nsamp` is the number of samples per event.
    vecs_orientation : list
        orbital angular momenta; list of tuples with entries of dimension 
        `(Nsamp, 3)`, where `Nsamp` is the number of samples per event.
    vecs_total_location : list
        total angular momenta; list of tuples with entries of dimension 
        `(Nsamp, 3)`, where `Nsamp` is the number of samples per event.
    """
    vecs_location = []
    vecs_orientation = []
    vecs_total_orientation = []
    for e, samples in tqdm(rwt_samples_dict.items()):
        df_e = pd.DataFrame(samples[ks1])
        df_e['sindec'] = np.sin(df_e['dec'])
        df_e['cosiota'] = np.cos(df_e['iota'])
        # assume reference frequency was default if not given
        if 'fref' not in df_e.columns:
            df_e['fref'] = fref
        # the following will apply the desired functions to each row
        # of the DataFrame
        nhats = df_e[ks2].apply(lambda x: get_location_vector(*x), axis=1)
        lhats = df_e[ks2].apply(lambda x: get_orientation_vector(*x), axis=1)
        jhats = df_e[ks3].apply(lambda x: get_total_orientation_vector(*x), axis=1)
        vecs_location.append(np.stack(nhats.values))
        vecs_orientation.append(np.stack(lhats.values))
        vecs_total_orientation.append(np.stack(jhats.values))
    return vecs_location, vecs_orientation, vecs_total_orientation

def get_sel_vectors(df_sel, fref=20, rng=np.random.default_rng()):
    """Get location (N) and orientation (J, L) vectors in Cartesian celestial 
    coordinates for a number of detected injections to evaluate the selection
    function.
    
    Arguments
    ---------
    df_sel : pd.DataFrame
        DataFrame containing parameters of detected injections as columns.
    fref : float
        reference frequency for spin orientations.
    rng : np.random.Generator
        random number generator used to generate random `phase` if that 
        parameter is not provided (opt.)
        
    Returns
    -------
    sel_vecs_n : pd.DataFrame
        binary locations; shape `(Ninj, 3)` for `Ninj` injections.
    sel_vecs_l : pd.DataFrame
        orbital angular momenta; shape `(Ninj, 3)` for `Ninj` injections.
    sel_vecs_j : pd.DataFrame
        total angu momenta; shape `(Ninj, 3)` for `Ninj` injections.
    """
    df_sel['mass_1'] = df_sel['m1'].values
    df_sel['mass_2'] = (df_sel['m1']*df_sel['q']).values
    if 'fref' not in df_sel.columns:
        df_sel['fref'] = fref
    if 'phase' not in df_sel.columns:
        df_sel['phase'] = rng.uniform(0, 2*np.pi, len(df_sel))
    sel_vecs_n = df_sel[ks2].apply(lambda x: get_location_vector(*x), axis=1)
    sel_vecs_l = df_sel[ks2].apply(lambda x: get_orientation_vector(*x), axis=1)
    sel_vecs_j = df_sel[ks3].apply(lambda x: get_total_orientation_vector(*x), axis=1)
    return sel_vecs_n, sel_vecs_l, sel_vecs_j

