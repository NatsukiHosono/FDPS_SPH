# This is a python script that converts u(rho, T), P(rho, T), Cs(rho,T), S(rho, T)
# to T(rho, u), P(rho, u), Cs(rho, u), S(rho, u), which is more useful for SPH calculaitons
#

import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.interpolate import interp1d
from scipy import interpolate

#----- A user has to change these three parameters  ----------------

inputfilename="granite.table.txt"    #input ANEOS file. This follows the format from iSALE
outputfilename="granite.rho_u.txt" #output ANEOS file
nu=120 #number of the grid for the internal energy (exponential)

#-------------------------------------------------------------------

# This function is to correct the original ANEOS format that does not include "E"
# This seems to  occur when the exponent reaches -101
def reformat(number):
    if number.find('E') == -1:
        exponent = "-101"
        mantissa  = number.split(exponent)      
        return float(mantissa[0])*10**float(exponent)
    else:
        mantissa, exponent= number.split('E')

    return float(mantissa)*10**float(exponent)


aneosfile= [line.split() for line in open(inputfilename)]

temperature=np.zeros(shape=(0,0))
density=np.zeros(shape=(0,0))

for i in range(1,len(aneosfile)):
    try:
        temperature=np.append(temperature,reformat(aneosfile[i][1]))
    except IndexError:
        nt=i-1
        break
    
for i in range(1,len(aneosfile), nt+1):
    density=np.append(density,reformat(aneosfile[i][0]))


nr=len(density) #density grid number

energy=np.zeros(shape=(nr,nt)) #J/kg
pressure=np.zeros(shape=(nr,nt)) #Pa
soundspeed=np.zeros(shape=(nr,nt)) #m/s
entropy=np.zeros(shape=(nr,nt)) #J/kg/K

i=1
for m in range(0,nr):
    for n in range(0,nt):
        try:
            energy[m][n]=reformat(aneosfile[i][2])
            pressure[m][n]=reformat(aneosfile[i][3])
            soundspeed[m][n]=reformat(aneosfile[i][4])
            entropy[m][n]=reformat(aneosfile[i][5])
            
        except IndexError: #skipping a line
            i=i+1
            energy[m][n]=reformat(aneosfile[i][2])
            pressure[m][n]=reformat(aneosfile[i][3])
            soundspeed[m][n]=reformat(aneosfile[i][4])
            entropy[m][n]=reformat(aneosfile[i][5])
        i=i+1


# Taking the min and max internal energy from the original ANEOS data
umin=np.min(energy)
umax=np.max(energy)

delta=(umax/umin)**(1.0/(nu-1))

new_energy=np.zeros(shape=(0,0))
for m in range(0,nu):
    new_energy=np.append(new_energy,umin*delta**m) #exponential grid

new_temperature=np.zeros(shape=(nr,nu))
new_pressure=np.zeros(shape=(nr,nu))
new_soundspeed=np.zeros(shape=(nr,nu))
new_entropy=np.zeros(shape=(nr,nu))


# 1D interpolation & extrapolation (linear)
for m in range(0,nu):

    # internal energy
    f_temperature = interpolate.interp1d(energy[m,:], temperature, kind='linear', fill_value='extrapolate')
    new_temperature[m][:]=f_temperature(new_energy)

    # pressure
    f_pressure  = interpolate.interp1d(temperature, pressure[m,:], kind='linear', fill_value='extrapolate')
    new_pressure[m][:]=f_pressure(new_temperature[m][:])

    # sound speed
    f_soundspeed  = interpolate.interp1d(temperature, soundspeed[m,:], kind='linear', fill_value='extrapolate')
    new_soundspeed[m][:]=f_soundspeed(new_temperature[m][:])

    # entropy
    f_entropy  = interpolate.interp1d(temperature, entropy[m,:], kind='linear', fill_value='extrapolate')
    new_entropy[m][:]=f_entropy(new_temperature[m][:])    


# producing a few output images to make sure that this fitting is doing an okay job
for m in range(0,nr, int(nr/6)):

    ax=[0, 0, 0, 0]

    fig=plt.figure(figsize=(10,6.128))

    ax[0]=fig.add_subplot(221)
    ax[1]=fig.add_subplot(222)
    ax[2]=fig.add_subplot(223)
    ax[3]=fig.add_subplot(224)

    ax[0].semilogy(temperature*1e-3, energy[m,:]*1e-6, '--', label = "original ANEOS")
    ax[0].semilogy(new_temperature[m, :]*1e-3, new_energy[:]*1e-6, '-.', label="modified")
    ax[1].semilogy(temperature*1e-3, pressure[m,:]*1e-6,'--', new_temperature[m, :]*1e-3, new_pressure[m, :]*1e-6,'-.')
    ax[2].plot(temperature*1e-3, soundspeed[m,:]*1e-3,'--', new_temperature[m, :]*1e-3, new_soundspeed[m, :]*1e-3,'-.')
    ax[3].plot(temperature*1e-3, entropy[m,:]*1e-3,'--', new_temperature[m, :]*1e-3, new_entropy[m, :]*1e-3,'-.')

    ax[0].legend(frameon=False)

    ax[0].set_ylabel('Energy (MJ/kg)',fontsize=10); ax[1].set_ylabel('Pressure (MPa)',fontsize=10); ax[2].set_ylabel('Sound Speed (km/s)',fontsize=10)
    ax[3].set_ylabel('Entropy (kJ/K/kg)',fontsize=10)
    ax[2].set_xlabel('Temperature ($10^3$ K)',fontsize=10); ax[3].set_xlabel('Temperature ($10^3$ K)',fontsize=10)
    
    fig.suptitle("Density: %3.3f kg/m$^3$" %(density[m]))
    fig.savefig("Density" + str(m) + ".png")

h = open(outputfilename,'w')
h.write("%i %i %s \n" % (nr, nu , ':Grid numbers for density and internal energy'))
h.write('Density (km/m3), Internal energy (kJ/kg), Temperature (K), Sound speed (m/s), Entropy (J/K/kg) \n')


for m in range(0,nr):
    for n in range(0,nu):
        h.write("%15.8E %15.8E %15.8E %15.8E %15.8E \n" % (density[m], new_energy[n], new_temperature[m][n],  new_soundspeed[m][n], new_entropy[m][n]))
