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
outputfilename="granite.rho_u.general.txt" #output ANEOS file
nu_new=120 #number of grid for new internal energy (exponential)
nr_new=240 #number of grid for new density

#-------------------------------------------------------------------

# This function is to correct the original ANEOS format that does not include "E"
# This seems to  occur when the exponent reaches -101


# --- reading ANEOS data
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
        
#--- defining new density grid -------------

# Taking the min and max internal energy from the original ANEOS data
rmin=np.min(density)
rmax=np.max(density)

delta=(rmax/rmin)**(1.0/(nr_new-1))

density_new=np.zeros(shape=(0,0))
for m in range(0,nr_new):
    density_new=np.append(density_new,rmin*delta**m) #exponential grid

#--- defining new internal energy grid  -------------

umin=np.min(energy)
umax=np.max(energy)

delta=(umax/umin)**(1.0/(nu_new-1))

energy_new=np.zeros(shape=(0,0))
for m in range(0,nu_new):
    energy_new=np.append(energy_new,umin*delta**m) #exponential grid

temperature_new=np.zeros(shape=(nr_new,nu_new))
pressure_new=np.zeros(shape=(nr_new,nu_new))
soundspeed_new=np.zeros(shape=(nr_new,nu_new))
entropy_new=np.zeros(shape=(nr_new,nu_new))
energy_interpolate=np.zeros(nt)

f_energy  = interpolate.interp2d(temperature, density, energy,kind='linear')
f_pressure  = interpolate.interp2d(temperature, density, pressure,kind='linear')
f_soundspeed  = interpolate.interp2d(temperature, density, soundspeed,kind='linear')
f_entropy  = interpolate.interp2d(temperature, density, entropy,kind='linear')


for m in range(0,nr_new):
    for n in range(0,nu_new):
        # at a fixed rho, this makes a relationship between energy and temperature
        energy_interpolate = f_energy(temperature, density_new[m])

        # temperature
        f_temperature = interpolate.interp1d(energy_interpolate, temperature, kind='linear', fill_value='extrapolate')
        temperature_new[m][n]=f_temperature(energy_new[n])

        # pressure
        pressure_new[m][n]=f_pressure(temperature_new[m][n],density_new[m])

        # sound speed
        soundspeed_new[m][n]=f_soundspeed(temperature_new[m][n],density_new[m])

        # entropy
        entropy_new[m][n]=f_entropy(temperature_new[m][n],density_new[m])



# producing a few output images to make sure that this fitting is doing an okay job
for m in range(0,nr, int(nr/6)):

    ax=[0, 0, 0, 0]

    fig=plt.figure(figsize=(10,6.128))

    ax[0]=fig.add_subplot(221)
    ax[1]=fig.add_subplot(222)
    ax[2]=fig.add_subplot(223)
    ax[3]=fig.add_subplot(224)

    for i in range(0,nr_new):
        if (density_new[i] >= density[m]):
            n=i
            break;

    
    ax[0].semilogy(temperature*1e-3, energy[m][:]*1e-6, '--', label = "original ANEOS")
    ax[0].semilogy(temperature_new[n][:]*1e-3, energy_new*1e-6, '-.', label="modified")

    
    ax[1].semilogy(temperature*1e-3, pressure[m,:]*1e-6,'--', temperature_new[n, :]*1e-3, pressure_new[n, :]*1e-6,'-.')
    ax[2].plot(temperature*1e-3, soundspeed[m,:]*1e-3,'--', temperature_new[n, :]*1e-3, soundspeed_new[n, :]*1e-3,'-.')
    ax[3].plot(temperature*1e-3, entropy[m,:]*1e-3,'--', temperature_new[n, :]*1e-3, entropy_new[n, :]*1e-3,'-.')

    ax[0].legend(frameon=False)

    ax[0].set_ylabel('Energy (MJ/kg)',fontsize=10); ax[1].set_ylabel('Pressure (MPa)',fontsize=10); ax[2].set_ylabel('Sound Speed (km/s)',fontsize=10)
    ax[3].set_ylabel('Entropy (kJ/K/kg)',fontsize=10)
    ax[2].set_xlabel('Temperature ($10^3$ K)',fontsize=10); ax[3].set_xlabel('Temperature ($10^3$ K)',fontsize=10)
    
    fig.suptitle("Density: %3.3f kg/m$^3$" %(density[m]))
    fig.savefig("Density" + str(m) + ".general.png")



h = open(outputfilename,'w')
h.write("%i %i %s \n" % (nr_new, nu_new, ':Grid numbers for density and internal energy'))
h.write('Density (km/m3), Internal energy (kJ/kg),Temperature (K),  Pressure (Pa), Sound speed (m/s), Entropy (J/K/kg) \n')


for m in range(0,nr_new):
    for n in range(0,nu_new):
        h.write("%15.8E %15.8E %15.8E %15.8E %15.8E %15.8E \n" % (density_new[m], energy_new[n], temperature_new[m][n],  pressure_new[m][n], soundspeed_new[m][n], entropy_new[m][n]))
