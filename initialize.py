import settings
import random
import math
import numpy as np
from numba import njit, prange
 
 
def InitializeAtoms():
    
    nx = 0
    ny = 0
    nz = 0
    n = 0
    x = np.zeros(shape=(settings.N))
    y = np.zeros(shape=(settings.N))
    z = np.zeros(shape=(settings.N))
    vx = np.zeros(shape=(settings.N))
    vy = np.zeros(shape=(settings.N))
    vz = np.zeros(shape=(settings.N))

    while nx < settings.n1:
        ny=0
        while ny < settings.n2:
            nz=0
            while nz < settings.n3*2:  # here we want to partilces
                if nz%2 == 0:
                    x0 = (nx+1/2) * settings.deltaxyz
                    y0 = (ny+1/2) * settings.deltaxyz
                    z0 = (nz+1/2-settings.bond_len/2) * settings.deltaxyz
                    
                    vx0 = 0.5 - random.randint(0, 1)
                    vy0 = 0.5 - random.randint(0, 1)
                    vz0 = 0.5 - random.randint(0, 1)
                        
                        
                    x[n] = x0
                    y[n] = y0
                    z[n] = z0
                    
                    vx[n] = vx0
                    vy[n] = vy0
                    vz[n] = vz0
                    n += 1
                        
                    nz += 1
                else:

                    x0 = (nx+1/2) * settings.deltaxyz
                    y0 = (ny+1/2) * settings.deltaxyz
                    z0 = (nz+1/2+settings.bond_len/2) * settings.deltaxyz
                    
                    vx0 = 0.5 - random.randint(0, 1)
                    vy0 = 0.5 - random.randint(0, 1)
                    vz0 = 0.5 - random.randint(0, 1)
                        
                        
                    x[n] = x0
                    y[n] = y0
                    z[n] = z0
                    
                    vx[n] = vx0
                    vy[n] = vy0
                    vz[n] = vz0
                    n += 1
                        
                    nz += 1
            ny +=1
        nx += 1
        
    # cancel the linear momentum
    svx = np.sum(vx)
    svy = np.sum(vy)
    svz = np.sum(vz)
    
    vx -= svx / settings.N 
    vy -= svy / settings.N 
    vz -= svz / settings.N 
    # svx = np.sum(vx)
    
    # rescale the velocity to the desired temperature
    Trandom = temperature(vx, vy, vz)
    vx, vy, vz = rescalevelocity(vx, vy, vz, settings.Tdesired, Trandom)
    
    # cancel the linear momentum
    svx = np.sum(vx)
    svy = np.sum(vy)
    svz = np.sum(vz)
    
    vx -= svx / settings.N 
    vy -= svy / settings.N 
    vz -= svz / settings.N 
    
    return x, y, z, vx, vy, vz

def temperature(vx, vy, vz):
    # receives units of [v] = nm/fs --> [v^2] = nm^2 /fs^2
    vsq = np.sum(vx*vx + vy*vy + vz*vz)

    # convert to kcal/mol:
    #   g/mol·(nm/fs)² → J/mol  by factor 1e9
    #   J/mol      → kcal/mol by dividing 4184
    conv = 1e9 / 4184.0

    K_kcal_per_mol = 0.5 * settings.mass * vsq * conv

    # equipartition: K = 3/2 N kB T  ⇒  T = 2K/(3N kB)
    return 2.0 * K_kcal_per_mol / (3.0 * settings.N * settings.kb)
    
def rescalevelocity(vx, vy, vz, T1, T2):
    fac = math.sqrt(T1 / T2)  # T1 is desired temperature
    vx = vx * fac
    vy = vy * fac
    vz = vz * fac
    return vx, vy, vz      

@njit(parallel=True)
def berendsen_thermostat(vx, vy, vz, T, T0, tau):
    """
    tau:    coupling strength
    T:      current Temperature
    T0:     Temp. to jump to / approach
    """
    lambdaa = np.sqrt(1 + settings.deltat / tau * (T0/T - 1))
    for i in prange(len(vx)):
        vx[i] *= lambdaa
        vy[i] *= lambdaa
        vz[i] *= lambdaa

    # return vx, vy, vz
    # no need to return anything since we just update the velocities


@njit(parallel=True)
def andersen_thermostat(vx, vy, vz, T0, nu):
    # vx, vy, vz : arrays of length N (particle velocities)
    # T0          : target temperature
    # nu           : collision frequency (collisions per unit time)

    p_collision = nu * settings.deltat   # probability each particle collides in this step
    # Precompute the Maxwell–Boltzmann velocity scale
    std = np.sqrt(settings.kb * T0 / settings.mass)

    for i in prange(settings.N):
        rdm = np.random.rand()
        if rdm < p_collision:
            # single components are normal distributed
            vx[i] = std * np.random.randn()
            vy[i] = std * np.random.randn()
            vz[i] = std * np.random.randn()

    # no need to return anything since we just update the velocities

if __name__ == '__main__':
    settings.init()
    #print(InitializeAtoms()[0])
    print(settings.deltaxyz)
    print(settings.bond_len)
    import matplotlib.pyplot as plt
    x, y, z, _, _, _ = InitializeAtoms()
    # plt.figure(figsize=[10,20])
    plt.scatter(y,z)
    plt.xlim([0, settings.l])
    plt.ylim([0, settings.l*2])
    plt.show()                  
    
    
    
