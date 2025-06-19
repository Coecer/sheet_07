import numpy as np
import settings
import initialize
import update
import force
import misc
from tqdm import tqdm
from numba import njit, prange

def SimpleSim(name, everyN):
    settings.init()
    x, y, z, vx, vy, vz = initialize.InitializeAtoms()
    vx = np.zeros(np.shape(vy))
    vy = np.zeros(np.shape(vy))
    vz = np.zeros(np.shape(vy))
    rlistlist= np.zeros((settings.eqsteps+1, settings.N//2))
    fx, fy, fz, rlistlist[0] = force.forceLJ(x, y, z, settings.xlo, settings.xhi, settings.ylo, settings.yhi, settings.zlo,settings.zhi,
                                     settings.eps, settings.sig, settings.cutoff, settings.bond_strength, settings.bond_len)
    

    fileoutputeq = open(name+ str(everyN)+ '_berendsen_eq', "w")
    misc.WriteTrajectory3d(fileoutputeq, 0,x,y,z)

    for i in tqdm(range(settings.eqsteps)):
        x, y, z, vx, vy, vz, fx, fy, fz, rlistlist[i+1] = update.VelocityVerlet(x, y, z, vx, vy, vz, fx, fy, fz, settings.xlo, settings.xhi, settings.ylo, settings.yhi,
                                                                    settings.zlo, settings.zhi, settings.eps, settings.sig, 
                                                                    settings.cutoff, settings.deltat, settings.mass, settings.bond_strength, settings.bond_len)


        # save shit every n
        if i % everyN == 0:

            misc.WriteTrajectory3d(fileoutputeq, 0,x,y,z)
    return rlistlist














































# def SimulationSheet7a(wT,Tname, thermostat_type='_berendsen'):
    
#     x, y, z, = np.array(0, 0), np.array(0, 0), np.array(-settings.bond_len_t0/2, settings.bond_len_t0/2) 
#     vx, vy, vz = np.array(0, 0), np.array(0, 0), np.array(0, 0)

#     fx, fy, fz, epot = force.forceLJ_and_Vibration(x, y, z, settings.xlo, settings.xhi, 
#                                                            settings.ylo, settings.yhi, settings.zlo,settings.zhi,
#                                                            settings.eps, settings.sig, settings.cutoff, settings.beta)
    

# def SimulationSheet7(wT,Tname, everyN, nu, rho, thermostat_type='_berendsen'):
#     # initialize
#     # print(f"Before:\n rho = {settings.rho:.2e}")
#     # print(f"L = {settings.l}")

#     if rho != settings.rho:
#         settings.rho = rho
#         settings.l = (settings.N/rho)**(1/3)
#         settings.xhi = settings.l
#         settings.yhi = settings.l
#         settings.zhi = settings.l        
#         settings.deltaxy = (settings.xhi - settings.xlo) / settings.n1
#         settings.deltaz  = (settings.zhi - settings.zlo) / settings.n3

#     # print(f"After:\n rho = {settings.rho:.2e}")
#     # print(f"L = {settings.l}")


#     x, y, z, vx, vy, vz = initialize.InitializeAtoms()
#     fx, fy, fz, epot, virial = force.forceLJ(x, y, z, settings.xlo, settings.xhi, settings.ylo, settings.yhi, settings.zlo,settings.zhi,
#                                      settings.eps, settings.sig, settings.cutoff, settings.beta)
    

    
#     print(f'vx, vy, vz = {vx[0], vy[0], vz[0]}')
#     print(f'fx, fy, fz = {fx[0], fy[0], fz[0]}')

#     if thermostat_type=='berendsen':

#         # open documents for eq run
#         if wT:
#             fileoutputeq = open(Tname+ str(everyN)+ '_berendsen_eq', "w")
#             fileoutputEn = open(Tname+ str(everyN)+ '_berendsen_Energyeq', "w")
#             misc.WriteTrajectory3d(fileoutputeq, 0,x,y,z) # noch anpassen an L, Lz
#             misc.WriteEnergy(fileoutputEn, 0, epot, update.KineticEnergy(vx, vy, vz,settings.mass),
#                               np.sum(vx**2), np.sum(vy**2), np.sum(vz**2), virial, initialize.temperature(vx, vy, vz))

#         temperatures_prod = np.ones(settings.nsteps)   
#         temperatures_eq = np.ones(settings.eqsteps)
#         T_curr = initialize.temperature(vx, vy, vz)
#         vx, vy, vz = initialize.rescalevelocity(vx, vy, vz, settings.Tdesired, T_curr)

#         # do eq run
#         for i in tqdm(range(settings.eqsteps)):
#             x, y, z, vx, vy, vz, fx, fy, fz, epot = update.VelocityVerlet(x, y, z, vx, vy, vz, fx, fy, fz, settings.xlo, settings.xhi, settings.ylo, settings.yhi,
#                                                                         settings.zlo, settings.zhi, settings.eps, settings.sig, 
#                                                                         settings.cutoff, settings.deltat, settings.mass, settings.beta)

#             # Temp rescaling
#             T_curr = initialize.temperature(vx, vy, vz)
#             initialize.berendsen_thermostat(vx, vy, vz, T_curr, settings.Tdesired, settings.tau)
#             # T_after = initialize.temperature(vx, vy, vz)
#             # if i < 10:
#             #     print(f"[step {i}] T_before={T_curr:.1f} K, T_after={T_after:.1f} K, λ={(T_after/T_curr):.5f}")
#             # temperatures_eq[i] = T_after

#             # save shit every n
#             if i % everyN == 0:
#                 # print(f"Current temperature = {temperatures_prod[i]:.1f} of step [{i}/[{settings.nsteps}]]")
#                 if wT:
#                     misc.WriteTrajectory3d(fileoutputeq, 0,x,y,z) # noch anpassen an L, Lz
#                     misc.WriteEnergy(fileoutputEn, i, epot, update.KineticEnergy(vx, vy, vz,settings.mass), 
#                                      np.sum(vx**2), np.sum(vy**2), np.sum(vz**2), virial, initialize.temperature(vx, vy, vz))

#         #open shit for prod
#         if wT:
#             fileoutputprod = open(Tname+ str(everyN)+ '_berendsen_prod', "w")
#             fileoutputEnprod = open(Tname+ str(everyN)+ '_berendsen_Energyprod', "w")
#             misc.WriteTrajectory3d(fileoutputprod, 0,x,y,z) # noch anpassen an L, Lz
#             misc.WriteEnergy(fileoutputEnprod, 0, epot, update.KineticEnergy(vx, vy, vz,settings.mass), 
#                               np.sum(vx**2), np.sum(vy**2), np.sum(vz**2), virial, initialize.temperature(vx, vy, vz))

#         temperatures_prod = np.ones(settings.nsteps)    

#         # do prod run
#         Ngr = 0
#         nbins = int(settings.l/2 / settings.dr)
#         hist = np.zeros(nbins) 
#         xlo,xhi,ylo,yhi,zlo,zhi = settings.xlo,settings.xhi,settings.ylo,settings.yhi,settings.zlo,settings.zhi
#         for i in tqdm(range(settings.nsteps)):  # changed this here to start from 0 (Jonas)
#             x, y, z, vx, vy, vz, fx, fy, fz, epot = update.VelocityVerlet(x, y, z, vx, vy, vz, fx, fy, fz, 
#                                                                         settings.xlo, settings.xhi, settings.ylo, settings.yhi, settings.zlo, settings.zhi,
#                                                                             settings.eps, settings.sig, settings.cutoff, settings.deltat, settings.mass, settings.beta)

#             # Temp rescaling
#             T_curr = initialize.temperature(vx, vy, vz)
#             # desired Temperature now --> Tjump

#             initialize.berendsen_thermostat(vx, vy, vz, T_curr, settings.Tdesired, settings.tau) 
#             # temperatures_prod[i] = initialize.temperature(vx, vy, vz)


#             # save shit every n
#             if i % everyN == 0:
#                 # print(f"Current temperature = {temperatures_prod[i]:.1f} of step [{i}/[{settings.nsteps}]]")
#                 if wT:
#                     misc.WriteTrajectory3d(fileoutputprod, 0,x,y,z) # noch anpassen an L, Lz
#                     misc.WriteEnergy(fileoutputEnprod, i, epot, update.KineticEnergy(vx, vy, vz,settings.mass),
#                                      np.sum(vx**2), np.sum(vy**2), np.sum(vz**2), virial, initialize.temperature(vx, vy, vz))
                    

#                 hist = update.update_hist(hist, x,y,z,
#                                           settings.xlo,settings.xhi,settings.ylo,settings.yhi,settings.zlo,settings.zhi,
#                                           settings.dr, settings.N, settings.l)
#                 Ngr += 1 # another position

#             g = update.calcg(Ngr,hist, settings.dr)

                    
#         return g
    


# def SimluationSheet5(wT,Tname, everyN, nu, thermostat_type='andersen'):
#     # initialize shit
#     settings.init()
#     x, y, z, vx, vy, vz = initialize.InitializeAtoms()
#     fx, fy, fz, epot, virial = force.forceLJ(x, y, z, settings.xlo, settings.xhi, settings.ylo, settings.yhi, settings.zlo,settings.zhi,
#                                      settings.eps, settings.sig, settings.cutoff)
    
#     if nu == settings.nu:
#         nu = settings.nu
#     else:
#         nu = nu


#     if thermostat_type=='berendsen':

#         # open documents for eq run
#         if wT:
#             fileoutputeq = open(Tname+ str(everyN)+ '_berendsen_eq', "w")
#             fileoutputEn = open(Tname+ str(everyN)+ '_berendsen_Energyeq', "w")
#             misc.WriteTrajectory3d(fileoutputeq, 0,x,y,z) # noch anpassen an L, Lz
#             misc.WriteEnergy(fileoutputEn, 0, epot, update.KineticEnergy(vx, vy, vz,settings.mass), np.sum(vx**2), np.sum(vy**2), np.sum(vz**2))

#         temperatures_prod = np.ones(settings.nsteps)   
#         temperatures_eq = np.ones(settings.eqsteps)
#         T_curr = initialize.temperature(vx, vy, vz)
#         vx, vy, vz = initialize.rescalevelocity(vx, vy, vz, settings.Tdesired, T_curr)

#         # do eq run
#         for i in tqdm(range(settings.eqsteps)):
#             x, y, z, vx, vy, vz, fx, fy, fz, epot = update.VelocityVerlet(x, y, z, vx, vy, vz, fx, fy, fz, settings.xlo, settings.xhi, settings.ylo, settings.yhi,
#                                                                         settings.zlo, settings.zhi, settings.eps, settings.sig, 
#                                                                         settings.cutoff, settings.deltat, settings.mass)
        
#             # Temp rescaling
#             T_curr = initialize.temperature(vx, vy, vz)
#             initialize.berendsen_thermostat(vx, vy, vz, T_curr, settings.Tdesired, settings.tau)
#             T_after = initialize.temperature(vx, vy, vz)
#             # if i < 10:
#             #     print(f"[step {i}] T_before={T_curr:.1f} K, T_after={T_after:.1f} K, λ={(T_after/T_curr):.5f}")
#             temperatures_eq[i] = T_after

#             # save shit every n
#             if i % everyN == 0:
#                 # print(f"Current temperature = {temperatures_prod[i]:.1f} of step [{i}/[{settings.nsteps}]]")
#                 if wT:
#                     misc.WriteTrajectory3d(fileoutputeq, 0,x,y,z) # noch anpassen an L, Lz
#                     misc.WriteEnergy(fileoutputEn, i, epot, update.KineticEnergy(vx, vy, vz,settings.mass), np.sum(vx**2), np.sum(vy**2), np.sum(vz**2))



#         #open shit for prod
#         if wT:
#             fileoutputprod = open(Tname+ str(everyN)+ '_berendsen_prod', "w")
#             fileoutputEnprod = open(Tname+ str(everyN)+ '_berendsen_Energyprod', "w")
#             misc.WriteTrajectory3d(fileoutputprod, 0,x,y,z) # noch anpassen an L, Lz
#             misc.WriteEnergy(fileoutputEnprod, 0, epot, update.KineticEnergy(vx, vy, vz,settings.mass), np.sum(vx**2), np.sum(vy**2), np.sum(vz**2))
        
#         temperatures_prod = np.ones(settings.nsteps)    

#         # do prod run
#         for i in tqdm(range(settings.nsteps)):  # changed this here to start from 0 (Jonas)
#             x, y, z, vx, vy, vz, fx, fy, fz, epot = update.VelocityVerlet(x, y, z, vx, vy, vz, fx, fy, fz, 
#                                                                         settings.xlo, settings.xhi, settings.ylo, settings.yhi, settings.zlo, settings.zhi,
#                                                                             settings.eps, settings.sig, settings.cutoff, settings.deltat, settings.mass)

#             # Temp rescaling
#             T_curr = initialize.temperature(vx, vy, vz)
#             # desired Temperature now --> Tjump
#             if i == settings.njump:
#                 settings.Tdesired = settings.Tjump

#             initialize.berendsen_thermostat(vx, vy, vz, T_curr, settings.Tdesired, settings.tau) 
#             temperatures_prod[i] = initialize.temperature(vx, vy, vz)


#             # save shit every n
#             if i % everyN == 0:
#                 # print(f"Current temperature = {temperatures_prod[i]:.1f} of step [{i}/[{settings.nsteps}]]")
#                 if wT:
#                     misc.WriteTrajectory3d(fileoutputeq, 0,x,y,z) # noch anpassen an L, Lz
#                     misc.WriteEnergy(fileoutputEn, i, epot, update.KineticEnergy(vx, vy, vz,settings.mass), np.sum(vx**2), np.sum(vy**2), np.sum(vz**2))

#         return temperatures_eq, temperatures_prod   
    

#     ########### b) do the same for the andersen thermostat ###########
#     elif thermostat_type=='andersen':
#         temperatures_eq = np.ones(settings.eqsteps)
#         temperatures_prod = np.ones(settings.nsteps)   

#         # rescale velocity once with  
#         T_curr = initialize.temperature(vx, vy, vz)
#         vx, vy, vz = initialize.rescalevelocity(vx, vy, vz, settings.Tdesired, T_curr) # get back to 300K

#         # open documents for eq run
#         if wT:
#             fileoutputeq = open(Tname+ str(everyN)+ '_andersen_eq', "w")
#             fileoutputEn = open(Tname+ str(everyN)+ '_andersen_Energyeq', "w")
#             misc.WriteTrajectory3d(fileoutputeq, 0,x,y,z) # noch anpassen an L, Lz
#             misc.WriteEnergy(fileoutputEn, 0, epot, update.KineticEnergy(vx, vy, vz,settings.mass), np.sum(vx**2), np.sum(vy**2), np.sum(vz**2))


#         # do eq run
#         for i in tqdm(range(settings.eqsteps)):
#             x, y, z, vx, vy, vz, fx, fy, fz, epot = update.VelocityVerlet(x, y, z, vx, vy, vz, fx, fy, fz, settings.xlo, settings.xhi, settings.ylo, settings.yhi,
#                                                                         settings.zlo, settings.zhi, settings.eps, settings.sig, 
#                                                                         settings.cutoff, settings.deltat, settings.mass)

#             # Temp rescaling
#             T_curr = initialize.temperature(vx, vy, vz)
#             initialize.andersen_thermostat(vx, vy, vz, settings.Tdesired, nu)
#             T_after = initialize.temperature(vx, vy, vz)
#             # if i < 10:
#             #     print(f"[step {i}] T_before={T_curr:.1f} K, T_after={T_after:.1f} K, λ={(T_after/T_curr):.5f}")
#             temperatures_eq[i] = T_after

#             if i % everyN:
#                 # print(f"Current temperature = {temperatures_prod[i]:.1f} of step [{i}/[{settings.nsteps}]]")
#                 if wT:
#                     misc.WriteTrajectory3d(fileoutputeq, 0,x,y,z) # noch anpassen an L, Lz
#                     misc.WriteEnergy(fileoutputEn, i, epot, update.KineticEnergy(vx, vy, vz,settings.mass), np.sum(vx**2), np.sum(vy**2), np.sum(vz**2))


#         #open documents for production run
#         if wT:
#             fileoutputprod = open(Tname+ str(everyN)+ '_andersen_prod', "w")
#             fileoutputEnprod = open(Tname+ str(everyN)+ '_andersen_Energyprod', "w")
#             misc.WriteTrajectory3d(fileoutputprod, 0,x,y,z) # noch anpassen an L, Lz
#             misc.WriteEnergy(fileoutputEnprod, 0, epot, update.KineticEnergy(vx, vy, vz,settings.mass), np.sum(vx**2), np.sum(vy**2), np.sum(vz**2))
        

#         # do prod run
#         for i in tqdm(range(settings.nsteps)):
#             x, y, z, vx, vy, vz, fx, fy, fz, epot = update.VelocityVerlet(x, y, z, vx, vy, vz, fx, fy, fz, settings.xlo, settings.xhi, settings.ylo, settings.yhi,
#                                                                         settings.zlo, settings.zhi, settings.eps, settings.sig, 
#                                                                         settings.cutoff, settings.deltat, settings.mass)
    

#             if i == settings.njump:
#                 settings.Tdesired = settings.Tjump

#             initialize.andersen_thermostat(vx, vy, vz, settings.Tdesired, nu)
#             temperatures_prod[i] = initialize.temperature(vx, vy, vz)


#             if i % everyN:
#                 # print(f"Current temperature = {temperatures_prod[i]:.1f} of step [{i}/[{settings.nsteps}]]")
#                 if wT:
#                     misc.WriteTrajectory3d(fileoutputprod, 0,x,y,z) # noch anpassen an L, Lz
#                     misc.WriteEnergy(fileoutputEnprod, i, epot, update.KineticEnergy(vx, vy, vz,settings.mass), np.sum(vx**2), np.sum(vy**2), np.sum(vz**2))

#         return temperatures_eq, temperatures_prod





# # ##################### Part c) & d) ##########################

# def SimluationSheet5_cd(wT,Tname, everyN, tau, thermostat_type='berendsen'):
#     # initialize shit
#     settings.init()
#     x, y, z, vx, vy, vz = initialize.InitializeAtoms()
#     fx, fy, fz, epot, virial = force.forceLJ(x, y, z, settings.xlo, settings.xhi, settings.ylo, settings.yhi, settings.zlo,settings.zhi,
#                                      settings.eps, settings.sig, settings.cutoff)
    
#     if tau == settings.tau:
#         tau = settings.tau
#     else:
#         tau = tau

#     print(tau)

#     if thermostat_type=='berendsen':

#         # open documents for eq run
#         if wT:
#             fileoutputeq = open(Tname+ str(everyN)+ '_berendsen_eq', "w")
#             fileoutputEn = open(Tname+ str(everyN)+ '_berendsen_Energyeq', "w")
#             misc.WriteTrajectory3d(fileoutputeq, 0,x,y,z) # noch anpassen an L, Lz
#             misc.WriteEnergy(fileoutputEn, 0, epot, update.KineticEnergy(vx, vy, vz,settings.mass), np.sum(vx**2), np.sum(vy**2), np.sum(vz**2))

#         temperatures_prod = np.ones(settings.nsteps)   
#         temperatures_eq = np.ones(settings.eqsteps)
#         T_curr = initialize.temperature(vx, vy, vz)
#         vx, vy, vz = initialize.rescalevelocity(vx, vy, vz, settings.Tdesired, T_curr)

#         # do eq run
#         for i in tqdm(range(settings.eqsteps)):
#             x, y, z, vx, vy, vz, fx, fy, fz, epot = update.VelocityVerlet(x, y, z, vx, vy, vz, fx, fy, fz, settings.xlo, settings.xhi, settings.ylo, settings.yhi,
#                                                                         settings.zlo, settings.zhi, settings.eps, settings.sig, 
#                                                                         settings.cutoff, settings.deltat, settings.mass)
        
#             # Temp rescaling
#             T_curr = initialize.temperature(vx, vy, vz)
#             initialize.berendsen_thermostat(vx, vy, vz, T_curr, settings.Tdesired, tau)
#             T_after = initialize.temperature(vx, vy, vz)
#             # if i < 10:
#             #     print(f"[step {i}] T_before={T_curr:.1f} K, T_after={T_after:.1f} K, λ={(T_after/T_curr):.5f}")
#             temperatures_eq[i] = T_after

#             # save shit every n
#             if i % everyN == 0:
#                 # print(f"Current temperature = {temperatures_prod[i]:.1f} of step [{i}/[{settings.nsteps}]]")
#                 if wT:
#                     misc.WriteTrajectory3d(fileoutputeq, 0,x,y,z) # noch anpassen an L, Lz
#                     misc.WriteEnergy(fileoutputEn, i, epot, update.KineticEnergy(vx, vy, vz,settings.mass), np.sum(vx**2), np.sum(vy**2), np.sum(vz**2))
            

#         #open shit for prod
#         if wT:
#             fileoutputprod = open(Tname+ str(everyN)+ '_berendsen_prod', "w")
#             fileoutputEnprod = open(Tname+ str(everyN)+ '_berendsen_Energyprod', "w")
#             misc.WriteTrajectory3d(fileoutputprod, 0,x,y,z) # noch anpassen an L, Lz
#             misc.WriteEnergy(fileoutputEnprod, 0, epot, update.KineticEnergy(vx, vy, vz,settings.mass), np.sum(vx**2), np.sum(vy**2), np.sum(vz**2))
        
#         temperatures_prod = np.ones(settings.nsteps)    

#         # do prod run
#         Ngr = 0
#         nbins = int(settings.l/2 / settings.dr)
#         hist = np.zeros(nbins) 
#         xlo,xhi,ylo,yhi,zlo,zhi = settings.xlo,settings.xhi,settings.ylo,settings.yhi,settings.zlo,settings.zhi

#         for i in tqdm(range(settings.nsteps)):  # changed this here to start from 0 (Jonas)
#             x, y, z, vx, vy, vz, fx, fy, fz, epot = update.VelocityVerlet(x, y, z, vx, vy, vz, fx, fy, fz, 
#                                                                         settings.xlo, settings.xhi, settings.ylo, settings.yhi, settings.zlo, settings.zhi,
#                                                                             settings.eps, settings.sig, settings.cutoff, settings.deltat, settings.mass)

#             # Temp rescaling
#             T_curr = initialize.temperature(vx, vy, vz)
#             initialize.berendsen_thermostat(vx, vy, vz, T_curr, settings.Tdesired, tau) 
#             temperatures_prod[i] = initialize.temperature(vx, vy, vz)


#             # save shit every n
#             if i % everyN == 0:
#                 # print(f"Current temperature = {temperatures_prod[i]:.1f} of step [{i}/[{settings.nsteps}]]")
#                 if wT:
#                     misc.WriteTrajectory3d(fileoutputprod, 0,x,y,z) # noch anpassen an L, Lz
#                     misc.WriteEnergy(fileoutputEnprod, i, epot, update.KineticEnergy(vx, vy, vz,settings.mass), np.sum(vx**2), np.sum(vy**2), np.sum(vz**2))
            


#             if i % settings.nAnalyze == 0:
#                 hist = update.update_hist(hist, x,y,z,
#                                           settings.xlo,settings.xhi,settings.ylo,settings.yhi,settings.zlo,settings.zhi,
#                                           settings.dr, settings.N, settings.l)
#                 Ngr += 1 # another position

#             g = update.calcg(Ngr,hist, settings.dr)


#         return temperatures_eq, temperatures_prod, g   
    



# # ##################### Part e) ##########################

# def SimluationSheet5_e(wT,Tname, everyN, nu, thermostat_type='andersen'):

#     # initialize shit
#     settings.init()
#     x, y, z, vx, vy, vz = initialize.InitializeAtoms()
#     fx, fy, fz, epot, virial = force.forceLJ(x, y, z, settings.xlo, settings.xhi, settings.ylo, settings.yhi, settings.zlo,settings.zhi,
#                                      settings.eps, settings.sig, settings.cutoff)

#     if thermostat_type=='andersen':
#         temperatures_eq = np.ones(settings.eqsteps)
#         temperatures_prod = np.ones(settings.nsteps)   

#         # rescale velocity once with  
#         T_curr = initialize.temperature(vx, vy, vz)
#         vx, vy, vz = initialize.rescalevelocity(vx, vy, vz, settings.Tdesired, T_curr) # get back to 300K

#         # open documents for eq run
#         if wT:
#             fileoutputeq = open(Tname+ str(everyN)+ '_andersen_eq', "w")
#             fileoutputEn = open(Tname+ str(everyN)+ '_andersen_Energyeq', "w")
#             misc.WriteTrajectory3d(fileoutputeq, 0,x,y,z) # noch anpassen an L, Lz
#             misc.WriteEnergy(fileoutputEn, 0, epot, update.KineticEnergy(vx, vy, vz,settings.mass), np.sum(vx**2), np.sum(vy**2), np.sum(vz**2))


#         # do eq run
#         for i in tqdm(range(settings.eqsteps)):
#             x, y, z, vx, vy, vz, fx, fy, fz, epot = update.VelocityVerlet(x, y, z, vx, vy, vz, fx, fy, fz, settings.xlo, settings.xhi, settings.ylo, settings.yhi,
#                                                                         settings.zlo, settings.zhi, settings.eps, settings.sig, 
#                                                                         settings.cutoff, settings.deltat, settings.mass)

#             # Temp rescaling
#             T_curr = initialize.temperature(vx, vy, vz)
#             initialize.andersen_thermostat(vx, vy, vz, settings.Tdesired, nu)
#             T_after = initialize.temperature(vx, vy, vz)
#             # if i < 10:
#             #     print(f"[step {i}] T_before={T_curr:.1f} K, T_after={T_after:.1f} K, λ={(T_after/T_curr):.5f}")

#             if i % everyN:
#                 # print(f"Current temperature = {temperatures_prod[i]:.1f} of step [{i}/[{settings.nsteps}]]")
#                 if wT:
#                     misc.WriteTrajectory3d(fileoutputeq, 0,x,y,z) # noch anpassen an L, Lz
#                     misc.WriteEnergy(fileoutputEn, i, epot, update.KineticEnergy(vx, vy, vz,settings.mass), np.sum(vx**2), np.sum(vy**2), np.sum(vz**2))


#         #open documents for production run
#         if wT:
#             fileoutputprod = open(Tname+ str(everyN)+ '_andersen_prod', "w")
#             fileoutputEnprod = open(Tname+ str(everyN)+ '_andersen_Energyprod', "w")
#             misc.WriteTrajectory3d(fileoutputprod, 0,x,y,z) # noch anpassen an L, Lz
#             misc.WriteEnergy(fileoutputEnprod, 0, epot, update.KineticEnergy(vx, vy, vz,settings.mass), np.sum(vx**2), np.sum(vy**2), np.sum(vz**2))
        

#         # do prod run
#         Ngr = 0
#         nbins = int(settings.l/2 / settings.dr)
#         hist = np.zeros(nbins) 
#         xlo,xhi,ylo,yhi,zlo,zhi = settings.xlo,settings.xhi,settings.ylo,settings.yhi,settings.zlo,settings.zhi

#         for i in tqdm(range(settings.nsteps)):
#             x, y, z, vx, vy, vz, fx, fy, fz, epot = update.VelocityVerlet(x, y, z, vx, vy, vz, fx, fy, fz, settings.xlo, settings.xhi, settings.ylo, settings.yhi,
#                                                                         settings.zlo, settings.zhi, settings.eps, settings.sig, 
#                                                                         settings.cutoff, settings.deltat, settings.mass)

#             initialize.andersen_thermostat(vx, vy, vz, settings.Tdesired, nu)
#             temperatures_prod[i] = initialize.temperature(vx, vy, vz)


#             if i % everyN:
#                 # print(f"Current temperature = {temperatures_prod[i]:.1f} of step [{i}/[{settings.nsteps}]]")
#                 if wT:
#                     misc.WriteTrajectory3d(fileoutputprod, 0,x,y,z) # noch anpassen an L, Lz
#                     misc.WriteEnergy(fileoutputEnprod, i, epot, update.KineticEnergy(vx, vy, vz,settings.mass), np.sum(vx**2), np.sum(vy**2), np.sum(vz**2))

#             if i % settings.nAnalyze == 0:
#                 hist = update.update_hist(hist, x,y,z,
#                                           settings.xlo,settings.xhi,settings.ylo,settings.yhi,settings.zlo,settings.zhi,
#                                           settings.dr, settings.N, settings.l)
#                 Ngr += 1 # another position

#             g = update.calcg(Ngr,hist, settings.dr)


#         return temperatures_eq, temperatures_prod, g   

