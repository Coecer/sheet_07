import settings
import numpy as np
from numba import njit, prange


#@njit


# @njit(parallel=True)
# def compute_B2_LJ(eps, sig, beta, r_max=1000, dr=0.001):
#     """Compute B2 for a Lennard-Jones potential by direct integration."""
#     r = np.arange(0.5*dr, r_max, dr)  # Equally spaced radial points
#     # r = r[r > 0]  # Avoid division by zero
#     B2 = 0.0
    
#     for i in prange(len(r)):
#         sf2 = (sig / r[i]) ** 2
#         sf6 = sf2 ** 3
#         U = 4 * eps * sf6 * (sf6 - 1)  # LJ potential
#         integrand = (np.exp(-beta * U) - 1) * r[i] ** 2
#         B2 += integrand * dr
    
#     B2 *= -2 * np.pi  # Spherical integration factor
#     return B2

# def compute_Pressure_from_B2(rho, kB, T_desired, B2):
#     """Calculates the pressure from: P_theo = rho * kB * T[1 + B2 * rho].
#     Not taking into accoutn B3. """
    
#     return rho * kB * T_desired * (1.0 + B2 * rho)



# @njit(parallel=True)
# def forceLJ_and_Vibration(x, y, z, xlo, xhi, ylo, yhi, zlo, zhi, eps, sig, cutoff):
        
#     fx = np.zeros(shape=len(x))
#     fy = np.zeros(shape=len(x))
#     fz = np.zeros(shape=len(x))
#     N = len(x)
    
#     epot = 0

#     for i in range(N-1):
#         for j in range(i+1, N):
#             rijx = pbc(x[i], x[j], xlo, xhi)
#             rijy = pbc(y[i], y[j], ylo, yhi)
#             rijz = pbc(z[i], z[j], zlo, zhi)
            
#             r2 = rijx * rijx + rijy * rijy + rijz * rijz
#             r = np.sqrt(r2)
#             # Calculate B2 integrand: [exp(-βU(r)) - 1]
#             if r2 < cutoff * cutoff:
#                 sf2 = sig*sig / r2
#                 sf6 = sf2 * sf2 * sf2
#                 epot += 4.*eps*sf6*(sf6 - 1.)   # kcal/mol
                
#                 # ff = 48.*eps*sf6*(sf6 - 0.5)/r2 # kcal/mol/nm^2
#                 f_mag = 48*eps*sf6*(sf6 - 0.5) / r   # kcal/(mol·nm)

#                 # apply to components 
#                 fx[i] -= f_mag * (rijx/r) # kcal/nm/mol
#                 fy[i] -= f_mag * (rijy/r) # kcal/nm/mol
#                 fz[i] -= f_mag * (rijz/r) # kcal/nm/mol
#                 fx[j] += f_mag * (rijx/r) # kcal/nm/mol
#                 fy[j] += f_mag * (rijy/r) # kcal/nm/mol
#                 fz[j] += f_mag * (rijz/r) # kcal/nm/mol           

#     return fx, fy, fz, epot # units: kcal/mol/nm, kcal/mol, kcal/mol


#### unit of the force: (kcal/mole)/nm
@njit(parallel=True)
def forceLJ(x, y, z, xlo, xhi, ylo, yhi, zlo, zhi, eps, sig, cutoff, k, bo):
    
    fx = np.zeros(shape=len(x))
    fy = np.zeros(shape=len(x))
    fz = np.zeros(shape=len(x))
    N = len(x)
    
    epot = 0

    for i in range(N-1):
        for j in range(i+1, N):
            if (i%2 == 0 and j != i+1) or (i%2 ==1):
                rijx = pbc(x[i], x[j], xlo, xhi)
                rijy = pbc(y[i], y[j], ylo, yhi)
                rijz = pbc(z[i], z[j], zlo, zhi)
                
                r2 = rijx * rijx + rijy * rijy + rijz * rijz
                r = np.sqrt(r2)
                # Calculate B2 integrand: [exp(-βU(r)) - 1]
                if r2 < cutoff * cutoff:
                    sf2 = sig*sig / r2
                    sf6 = sf2 * sf2 * sf2
                    epot += 4.*eps*sf6*(sf6 - 1.)   # kcal/mol
                    
                    # ff = 48.*eps*sf6*(sf6 - 0.5)/r2 # kcal/mol/nm^2
                    f_mag = 48*eps*sf6*(sf6 - 0.5) / r   # kcal/(mol·nm)

                    # apply to components 
                    fx[i] -= f_mag * (rijx/r) # kcal/nm/mol
                    fy[i] -= f_mag * (rijy/r) # kcal/nm/mol
                    fz[i] -= f_mag * (rijz/r) # kcal/nm/mol
                    fx[j] += f_mag * (rijx/r) # kcal/nm/mol
                    fy[j] += f_mag * (rijy/r) # kcal/nm/mol
                    fz[j] += f_mag * (rijz/r) # kcal/nm/mol 
            else:
                rijx = pbc(x[i], x[j], xlo, xhi)
                rijy = pbc(y[i], y[j], ylo, yhi)
                rijz = pbc(z[i], z[j], zlo, zhi)

                r = np.sqrt(rijx**2 + rijy**2 + rijz**2)
                delta = r - bo
                f_mag = -k * delta / r  # normalized force direction

                f_x = f_mag * rijx
                f_y = f_mag * rijy
                f_z = f_mag * rijz

                fx[i] -= f_x
                fy[i] -= f_y
                fz[i] -= f_z

                fx[j] += f_x
                fy[j] += f_y
                fz[j] += f_z                






    
    return fx, fy, fz # units: kcal/mol/nm, kcal/mol, kcal/mol
  
@njit  
def pbc(xi, xj, xlo, xhi):
    
    l = xhi-xlo
    
    xi = xi % l
    xj = xj % l
    
    rij = xj - xi  
    if abs(rij) > 0.5*l:
        rij = rij - np.sign(rij) * l 
        
    return rij




# @njit
# def pbc(xi, xj, xlo, xhi):
#     l = xhi - xlo
#     rij = xi - xj
#     if rij > 0.5 * l:
#         rij -= l
#     elif rij < -0.5 * l:
#         rij += l
#     return rij
