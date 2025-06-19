#settings

# velocity: m*s^-1
# position: m
# acceleration: ms^-2
# energy: joule
# temperature: K

import numpy as np
from matplotlib import rcParams

# def reinitialize():  ## I assume we dont need this here?? -c
#     """"Updates all associated values when some
#       initialization values have changed"""
#     deltaxy = (xhi - xlo) / n1
#     deltaz  = (zhi - zlo) / n3
#     beta = 1 / (kb*Tdesired) 
#     rho = 0.05 * sig**(-3) # rho in units of nm^-3
#     l = (N/rho)**(1/3)
#     cutoff = 2.5*sig
#     N = n1*n2*n3
#     xhi = l
#     yhi = l
#     zhi = l
#     dr = sig/30
#     drxy = l/100

def init():
    global eqsteps
    eqsteps = 10         # != 50_000
    global nsteps            # number of time step to analyze
    nsteps = 10             # != 100_000
    # global njump                # number of time steps until temp. jump
    # njump = 400              # 40 000 actually 0 in prod run

    # global nAnalyze
    # nAnalyze = 50

    global mass              # mass of the LJ particles (gram/mol)
    mass = 13.93             # changed here to nitrogen mass

    global kb                # boltzmann's constant (kcal/mol/K) 
    kb = 0.0019858775
    # kb = 0.0019858775  * 4.184e-6
    global Tdesired          # temperature of the experiment in K
    Tdesired = 300.
    global beta                 # desired beat (kcal/mol)
    beta = 1 / (kb*Tdesired) 
    global tau
    tau = 100                # coupling parameter 
    # global nu
    # nu = 0.01             # collision frequency in (fs)^-1
    global eps               # eps in LJ (kcal/mol)
    eps = 0.297741315                  # (kcal/mol)
    # eps = 0.297741315 * 4.184e-6       # (kcal/mol)
    global sig              # r0 in LJ
    sig = 0.188             # in nm
    global cutoff            # cutoff arbitrary at 2.5 r0
    cutoff = 2.5*sig
    global deltat            # time step (fs)
    deltat = 2

    global bond_strength
    bond_strength = 9793        # kb/mol/nm^2
    global bond_len
    bond_len = 0.107            # b0 = 1.07 Angstroem
    global bond_len_t0
    bond_len_t0 = 0.11          # bond length at t=0
    # global mass_N
    # mass_N = 13.97              # mass of one nitrogen atom in g/mol


    
    # number of particle = n1*n2 distributed on s square lattice
    global n1
    n1 = 2
    global n2
    n2 = 2
    global n3
    n3 = 2
    global N
    N = n1*n2*n3*2

    # desired density
    global rho          # rho = N/V
    rho = 0.25 * sig**(-3) # rho in units of nm^-3 ## Achtung! ist moleküldichte nicht atomdichte!!
    global l
    l = (N/rho)**(1/3)

    # # for histogramms
    # global dr 
    # dr = sig/30
    # global drxy
    # drxy = l/100
    

    # box size
    global xlo
    xlo = 0.
    global xhi
    xhi = l
    global ylo
    ylo = 0.
    global yhi
    yhi = l
    global zlo
    zlo = 0.
    global zhi
    zhi = l
    
    global deltaxyz # lattice parameter to setup the initial configuration on a lattice
    deltaxyz = (xhi - xlo)/n1
    # global deltaz
    # deltaz = (zhi - zlo)/n3
    
    #rescaling of temperature
    global Trescale
    Trescale = 1 #1 = rescale temperature; 0 = no rescaling
    
    

    # Set global font family and sizes
    rcParams['font.family'] = 'Times New Roman'       # or 'Arial', 'Times New Roman', etc.
    rcParams['axes.labelsize'] = 18               # Axis label font size
    rcParams['xtick.labelsize'] = 13              # X-axis number size
    rcParams['ytick.labelsize'] = 13              # Y-axis number size
    rcParams['legend.fontsize'] = 14              # Legend font size
    rcParams['axes.titlesize'] = 20               # Title font size

    # Major tick size and width
    rcParams['xtick.major.size'] = 7
    rcParams['xtick.major.width'] = 1.5
    rcParams['ytick.major.size'] = 7
    rcParams['ytick.major.width'] = 1.5

    # Minor tick size and width
    rcParams['xtick.minor.size'] = 4
    rcParams['xtick.minor.width'] = 1
    rcParams['ytick.minor.size'] = 4
    rcParams['ytick.minor.width'] = 1

    # Optional: direction (in, out, or inout)
    rcParams['xtick.direction'] = 'in'
    rcParams['ytick.direction'] = 'in'

    # Disable scientific notation on all axes
    rcParams['axes.formatter.useoffset'] = False
    rcParams['axes.formatter.use_mathtext'] = False  # Optional: disables mathtext like 1e6 in fancy font


    # for plotting on white background, up to 14 colors
    color_list = [
        'tab:blue',
        'tab:orange',
        'tab:green',
        'tab:red',
        'tab:purple',
        'tab:brown',
        'tab:pink',
        'tab:gray',
        'tab:olive',
        'tab:cyan',
        'darkblue',    
        'mediumvioletred',
        'darkgreen'
    ]
    
