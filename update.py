# module containing integration sheme:
"""
     --------- units --------- 
     [r] = nm; [t] = fs; [epsilon] = kcal/mole; [m] = gram/mole;  [F] = (kcal/mole)/nm; 
     [T] = K; [v] = nm/fs; [k_B] = (kcal/mole)/K
     conversion factor:
         from kcal*fs*fs/gram/nm to nm: ???
         from kcal*fs/gram/nm to nm/fs: ??? 
"""

import settings
import force
import numpy as np
from numba import njit, prange

@njit(parallel=True)
def VelocityVerlet(x, y, z, vx, vy, vz, fx, fy, fz, xlo, xhi, ylo, yhi, zlo, zhi, eps, sigma, cutoff, deltat, mass, k, bo):

    # conversion factor
    fx0 = np.zeros(shape=len(x))
    fy0 = np.zeros(shape=len(x))
    fz0 = np.zeros(shape=len(x))
    N = len(x)
    dt = deltat
    #mass = mass
    
    #update the position at t+dt
    for i in prange(N):
        x[i] = (x[i] + vx[i] * dt + 4.184e-6 * fx[i] * dt * dt * 0.5 / mass ) % (xhi-xlo)
        y[i] = (y[i] + vy[i] * dt + 4.184e-6 * fy[i] * dt * dt * 0.5 / mass ) % (yhi-ylo)
        z[i] = (z[i] + vz[i] * dt + 4.184e-6 * fz[i] * dt * dt * 0.5 / mass ) % (zhi-zlo)
        
    # save the force at t
    fx0 = fx         # tested it also only with fx and it still worked
    fy0 = fy         # tested it also only with fy and it still worked
    fz0 = fz         # tested it also only with fz and it still worked
    # update acceleration at t+dt
    fx, fy, fz= force.forceLJ(x, y, z, xlo, xhi, ylo, yhi, zlo, zhi, eps, sigma, cutoff, k, bo)

    # conversion from [f_i/mass] = fs/nm * kcal/g --> 4.184*10^15 J/kg = 4.184*10^15 m/s^2 = 4.184e-6 nm/fs^2
    # conv_fac_ToMultiplyTo_vel = 4.184e-6               
    # update the velocity
    for i in prange(N):        
        vx[i] += 4.184e-6 * 0.5 * dt * (fx[i] + fx0[i]) / mass
        vy[i] += 4.184e-6 * 0.5 * dt * (fy[i] + fy0[i]) / mass
        vz[i] += 4.184e-6 * 0.5 * dt * (fz[i] + fz0[i]) / mass
    
    return x, y, z, vx, vy, vz, fx, fy, fz # units: nm, nm/fs, kcal/mol/nm, kcal/mol

 
@njit(parallel=True)
def KineticEnergy(vx, vy, vz, mass):
# calcualte the kinetic energy in joule
    ekin = 0
    N = len(vx)
    i = 0
    
    for i in prange(N):
        ekin += 0.5 * mass * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i]* vz[i])

    # ekin is in units of: [mass*vel^2] = g/mol * (nm/fs)^2 = 10^-3 * (10^-18/10^-30) kg*m^2/s^2 /mol 
    #  = 10^+9 J/mol = 10^+9 * 1/4184 = conv_fac_ToMultiplyTo_ekin
    conv_fac_ToMultiplyTo_ekin = 10**9 * 1/4184
    ekin = conv_fac_ToMultiplyTo_ekin * ekin
    return ekin

def compute_epot_from_gr(g_of_r, dr, rho, sigma, epsilon):
    """
    Given:
      - g_of_r:  1D array of length nbins, the time‐averaged RDF g(r) for r in [0, L/2)
      - dr:      bin width used to build g(r)
      - rho:     number density N/V
      - sigma:   LJ sigma
      - epsilon: LJ epsilon
    
    Returns:
      - E_pot_per_particle (float):  the ensemble‐average potential energy per particle 
        computed via 2πρ ∫ u_LJ(r) g(r) r^2 dr.
    """
    nbins = len(g_of_r)
    # bin centers r_i = (i + 0.5) * dr for i=0..nbins-1
    r_centers = (np.arange(nbins) + 0.5) * dr
    
    # LJ potential at r_i
    # (We assume r_centers never equals 0; dr>0 guarantees r_centers[0] = 0.5*dr > 0.)
    uLJ = 4.0 * epsilon * ( (sigma / r_centers)**12 - (sigma / r_centers)**6 )
    
    # integration
    E_pot = 2.0 * np.pi * rho * np.sum( uLJ * g_of_r * r_centers**2 * dr )
    return E_pot


######### added by Jonas ###########
def calcg(Ngr, hist, dr):
    r= r = np.arange(len(hist)) * dr    # hist goes up to l/2
    nid = 4*np.pi *settings.rho /3 * ((r+ dr)**3-r**3) # type: ignore
    n = hist/ settings.N /Ngr
    return n/nid


@njit(parallel=True)
def update_hist(hist, x,y,z,xlo,xhi,ylo,yhi,zlo,zhi, dr, N, l):
    l_half = (xhi - xlo) * 0.5
    for i in prange(N-1):
        for j in range(i+1,N):
            rijx = force.pbc(x[i], x[j], xlo, xhi) # calculate pbc distance
            rijy = force.pbc(y[i], y[j], ylo, yhi)
            rijz = force.pbc(z[i], z[j], zlo, zhi)
            r = np.sqrt(rijx**2 + rijy**2 + rijz**2)
            if r > 0.0 and r < l_half:
                bin = int(r / dr) # find the bin
                hist[bin] += 2 # we are counting pairs
    return hist


# def calc_rho(Ngr, hist, dr):
#     rho = hist / settings.N / Ngr / dr
#     return rho


# @njit
# def update_hists(hist_x, hist_y, hist_z, x,y,z, dr, drxy, N):
#     for i in prange(N):
#         binx = int(x[i]/drxy)
#         hist_x[binx] += 1
#         biny = int(y[i]/drxy)
#         hist_y[biny] += 1
#         binz = int(z[i]/dr)
#         hist_z[binz] += 1
#     return hist_x, hist_y, hist_z




    # for i in prange(N-1):
    #     # j = i + 1
    #     for j in range(i+1,N):
    #         xij = force.pbc(x[i], x[j], settings.xlo, settings.xhi) # calculate pbc distance
    #         yij = force.pbc(y[i], y[j], settings.ylo, settings.yhi) # can never exceed l/2
    #         zij = abs(z[i] - z[j])         # can never exceed 2l
 
    #         if xij != 0:
    #             bin = int(xij / dr) # find the bin
    #             hist_x[bin] += 2 # we are counting pairs
    #         if yij != 0:
    #             bin = int(yij / dr) # find the bin
    #             hist_y[bin] += 2 # we are counting pairs
    #         if zij != 0:
    #             bin = int(zij / dr) # find the bin
    #             hist_z[bin] += 2 # we are counting pairs

    # return hist_x, hist_y, hist_z