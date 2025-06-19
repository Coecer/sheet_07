import settings
import numpy as np
import math
import sys

def save_pressure_results(filename: str, pressure_means: np.ndarray,pressure_stds: np.ndarray,
                           pressure_std_standard_err_of_mean: float, P_theo: float) -> None:

    with open(filename, 'w') as f:
        f.write("# MD Simulation Pressure Analysis (Exercise 6.1 b (iii),(iv))\n")
        f.write(f"Theoretical pressure (from B2): P_theo = {P_theo:.6e}\n\n")
        for i, (mean, std) in enumerate(zip(pressure_means, pressure_stds), start=1):
            f.write(f"Block {i:1d}: <P{i}> = {mean:.6e} \u00b1 {std:.6e}\n")
        f.write(f"\nStandard error of the mean: {pressure_std_standard_err_of_mean:.6e}\n")
    
    print(f"Results written to {filename}")



def WriteTopology(filename, numb_atoms, atom_types, numb_bonds, x, y, z):
    """input has to be already have calculated bonds correctly
    """ 
    with open(filename, "w") as f:
        # block 1 -- number of atoms, atom types and bonds
        f.write("%i atoms\n" % numb_atoms)  # for atom_ids later
        f.write("%i atom types\n" % atom_types)     # in our case this should be 1 bc we only have N2
        f.write("%i bonds\n" % numb_bonds)          # numb_bonds != numb_atoms/2 in our case 
        f.write("1 bond types\n\n")      # always 1 for us here

        # Block 2 -- initial box coordinates
        # f.write("ITEM: BOX BOUNDS \n")
        f.write("%e %e xlo xhi \n" % (0, settings.l))
        f.write("%e %e xlo xhi \n" % (0, settings.l))
        f.write("%e %e xlo xhi \n\n" % (0, settings.l))

        # Block 3 -- Atom Block with atom_id, mol_id, atom_type
        numb_molecules = numb_atoms//numb_bonds  # rounds to lower value integer
        charge = 0.0
        atom_types_arr =  list(range(atom_types)) * numb_molecules
        f.write("Atoms\n\n") 
        atom_id = 0 
        charge = 0.0
        for mol_id in range(numb_bonds):  # each molecule = 1 bond = 2 atoms
            for _ in range(2):
                # third column = 1 always since we only have nitrogen N
                f.write("%i %i 1 %e %e %e %e\n" % (atom_id, mol_id,
                                                    charge, x[atom_id], y[atom_id], z[atom_id]))
                atom_id += 1
        # Block 4 -- Bonds
        f.write("\nBonds\n\n") 
        for bond_id, atom_id_even, atom_id_odd in zip(range(numb_bonds), range(0, numb_atoms, 2), range(1, numb_atoms, 2)):
            f.write("%i 1 %i %i\n" % (bond_id, atom_id_even, atom_id_odd))


def WriteEnergy(fileenergy, itime, epot, ekin, vx2, vy2, vz2, virial, Temp):
    
    fileenergy.write("%i %e %e %e %e %e %e %e\n" % (itime, epot, ekin, vx2, vy2, vz2, virial, Temp))

def WriteTrajectory(fileoutput, itime, x, y):

    fileoutput.write("ITEM: TIMESTEP \n")
    fileoutput.write("%i \n" % itime)
    fileoutput.write("ITEM: NUMBER OF ATOMS \n")
    fileoutput.write("%i \n" % (settings.n1*settings.n2))
    fileoutput.write("ITEM: BOX BOUNDS \n")
    fileoutput.write("%e %e xlo xhi \n" % (settings.xlo, settings.xhi))
    fileoutput.write("%e %e xlo xhi \n" % (settings.ylo, settings.yhi))
    fileoutput.write("%e %e xlo xhi \n" % (-1, 1))
    fileoutput.write("ITEM: ATOMS id type x y z \n")
    
    for i in range(0, len(x)):
        z = 0
        fileoutput.write("%i %i %e %e %e \n" % (i, i, (x[i]%settings.xhi), (y[i]%settings.yhi), z))
        

def inputset():
    return settings.xlo, settings.xhi, settings.ylo, settings.yhi, settings.eps, settings.sig, settings.cutoff, settings.deltat, settings.mass
    
def squarevelocity(vx, vy, mass):
    vx2 = 0
    vy2 = 0
    i = 0
    for i in range(0, len(vx)):
        vx2 += vx[i]**2
        vy2 += vy[i]**2
    return 0.5*mass*vx2, 0.5*mass*vy2


def WriteTrajectory3d(fileoutput, itime, x, y, z):

    fileoutput.write("ITEM: TIMESTEP \n")
    fileoutput.write("%i \n" % itime)
    fileoutput.write("ITEM: NUMBER OF ATOMS \n")
    fileoutput.write("%i \n" % (settings.N))
    fileoutput.write("ITEM: BOX BOUNDS \n")
    fileoutput.write("%e %e xlo xhi \n" % (0, settings.l))
    fileoutput.write("%e %e xlo xhi \n" % (0, settings.l))
    fileoutput.write("%e %e xlo xhi \n" % (0, settings.l))
    fileoutput.write("ITEM: ATOMS id type x y z \n")
    
    for i in range(0, settings.N):
        fileoutput.write("%i %i %e %e %e \n" % (i, i, x[i] % settings.l, y[i] % settings.l, z[i] % (settings.l)))
    
    
