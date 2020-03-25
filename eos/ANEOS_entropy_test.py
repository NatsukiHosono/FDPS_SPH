# This is a python script that converts u(rho, T), P(rho, T), Cs(rho,T), S(rho, T)
# to T(rho, u), P(rho, u), Cs(rho, u), S(rho, u), which is more useful for SPH calculaitons
#

import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.interpolate import interp1d, interp2d
from scipy import interpolate

#----- A user has to change these three parameters  ----------------
sphfilename = "results.granite.00100_00001_00000.dat"
inputfilename="granite.rho_u.txt"    #input ANEOS file. This follows the format from iSALE
outputfilename="output.txt"
particle_eos = 0 #0 for mantle 1 for core
#-------------------------------------------------------------------

aneosfile= [line.split() for line in open(inputfilename)]
nr = int(aneosfile[0][1])
nu = int(aneosfile[0][2])

temperature=np.zeros(shape=(nr,nu))
soundspeed=np.zeros(shape=(nr,nu))
pressure=np.zeros(shape=(nr,nu))
entropy=np.zeros(shape=(nr,nu))
density=np.zeros(shape=(nr))
energy=np.zeros(shape=(nu))

k=2
for m in range(0,nr):
    density[m]=aneosfile[k][0]
    for n in range(0,nu):
        energy[n]=aneosfile[k][1]
        temperature[m][n]=aneosfile[k][2]
        pressure[m][n]=aneosfile[k][3]
        soundspeed[m][n]=aneosfile[k][4]
        entropy[m][n]=aneosfile[k][5]
        k = k + 1

# making 2D interpolation.. thins transpose is weird but I think I have to do it
entropy_t = entropy.transpose()
pressure_t = pressure.transpose()
f_entropy = interpolate.interp2d(density, energy, entropy_t, kind='linear')
f_pressure = interpolate.interp2d(density, energy, pressure_t, kind='linear')


sphfile= [line.split() for line in open(sphfilename)]
n_SPH = len(sphfile)

density_SPH=np.zeros(shape=(n_SPH))
energy_SPH=np.zeros(shape=(n_SPH))
eos_SPH=np.zeros(shape=(n_SPH))
entropy_SPH=np.zeros(shape=(n_SPH))
rr_SPH=np.zeros(shape=(n_SPH))
PP_SPH=np.zeros(shape=(n_SPH))
pressure_SPH=np.zeros(shape=(n_SPH))
mass_SPH=np.zeros(shape=(n_SPH))

for m in range(2,n_SPH):
    eos_SPH[m]=int(sphfile[m][1])
    density_SPH[m] = float(sphfile[m][9])
    energy_SPH[m] = float(sphfile[m][10])
    rr_SPH[m] =  np.sqrt(float(sphfile[m][3])**2.0+ float(sphfile[m][4])**2.0+ float(sphfile[m][5])**2.0)*1e-6
    PP_SPH[m] =   float(sphfile[m][11])
    mass_SPH[m] =   float(sphfile[m][2])
    
    if eos_SPH[m]==particle_eos:
        entropy_SPH[m]=f_entropy(density_SPH[m], energy_SPH[m])
        pressure_SPH[m]=f_pressure(density_SPH[m], energy_SPH[m])        


h = open(outputfilename,'w')



h.write('# ID, radius, pressure (SPH, GPa), presure (interpolated, GPa), density (SPH), energy(SPH), entropy (SPH) \n')

for m in range(2,n_SPH):
    if eos_SPH[m]==particle_eos:
        h.write("%i  %15.8E  %15.8E   %15.8E %15.8E %15.8E  %15.8E \n" % (m, rr_SPH[m], PP_SPH[m]*1e-9, pressure_SPH[m]*1e-9, density_SPH[m], energy_SPH[m], entropy_SPH[m]))
      
        
