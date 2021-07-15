# This is a python script that converts u(rho, T), P(rho, T), Cs(rho,T), S(rho, T)
# to T(rho, u), P(rho, u), Cs(rho, u), S(rho, u), which is more useful for SPH calculaitons
#

import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.interpolate import interp1d, interp2d
from scipy import interpolate

#----- A user has to change these three parameters  ----------------

inputfilename="duniteS.table.txt"    #input ANEOS file. This follows the format from iSALE
outputfilename="duniteS2.rho_u.txt" #output ANEOS file
nu=120 #number of the grid for the internal energy (exponential)

#-------------------------------------------------------------------

# This function is to correct the original ANEOS format that does not include "E"
# This seems to  occur when the exponent reaches -101
def reformat(number):

    if number.find('E') == -1:       
        for i in range(100,200):
            number_to_find =  '-' + str(i)
            if number.find(number_to_find) >= 1:
                exponent = number_to_find


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

        #print(energy[m][n], pressure[m][n],soundspeed[m][n],entropy[m][n])
        

# Taking the min and max internal energy from the original ANEOS data
umin=np.min(energy)*1.01 # this is added to make interpolation working
umax=np.max(energy)
if umin < 0.0:
    umin = 1.0

delta=(umax/umin)**(1.0/(nu*1.0-1.0))

print(nu,nr,delta, umax,umin)
#sys.exit()


new_energy=np.zeros(shape=(0,0))
for m in range(0,nu):
    new_energy=np.append(new_energy,umin*delta**m) #exponential grid

new_temperature=np.zeros(shape=(nr,nu))
new_pressure=np.zeros(shape=(nr,nu))
new_soundspeed=np.zeros(shape=(nr,nu))
new_entropy=np.zeros(shape=(nr,nu))


def interpolation(p11, p12, p21, p22, alpha, beta):
    A1 = p11 + alpha*(p12-p11) #fixed density
    A2 = p21 + alpha*(p22-p21) #fixed density

    return A1 + beta*(A2-A1)

def rho_u(energy_list, density_list, temperature_list, u_in, rho_in):
    if rho_in < density[0]:
        print("too small density")
        sys.exit()

    if rho_in>=density[len(density)-1]:
        beta=0.0
        index0=len(density)-1
        index1=index0+1
        
    else:
        density_index = density_list.index(density[density>rho_in][0])
        beta = (rho_in - density_list[density_index-1])/(density_list[density_index]-density_list[density_index-1])

        index0=density_index-1
        index1=density_index+1
    

        
    temp_i_1 = [0, 0]
    n_i = [0, 0]
    
    k=0

    for i in range(index0, index1):
        energy_i_1 = energy[i][:]
        energy_list_i_1 = list(energy_i_1)

        if u_in >  energy_list_i_1[len(energy_i_1)-1]:
            temp_i_1[k]=temperature[len(energy_i_1)-1]
            
        else:
            n_i[k] = energy_list_i_1.index(energy_i_1[energy_i_1 >= u_in][0])

            if n_i[k]==0:
                alpha_i=0.0
                temp_i_1[k]=temperature[0]

            else:
                alpha_i = (u_in - energy_i_1[n_i[k]-1])/(energy_i_1[n_i[k]]-energy_i_1[n_i[k]-1])
                temp_i_1[k] = temperature[n_i[k]-1] + alpha_i * (temperature[n_i[k]] - temperature[n_i[k]-1])
        k=k+1

    if rho_in>=density[len(density)-1]:
        temp_out =  temp_i_1[0]
        n_den_s=len(density)-1
        n_den_l=len(density)-1
        
    else:
        temp_out = temp_i_1[0] +  beta * (temp_i_1[1]-temp_i_1[0])
        n_den_s=density_index-1
        n_den_l=density_index
        
    if  np.abs(temp_out-temperature[len(temperature)-1])<0.001:
        n_max =  len(energy_i_1)-1
        new_press= interpolation(pressure[n_den_s][n_max], pressure[n_den_s][n_max], pressure[n_den_l][n_max], pressure[n_den_l][n_max], 0.0, beta)
        new_entropy= interpolation(entropy[n_den_s][n_max], entropy[n_den_s][n_max], entropy[n_den_l][n_max], entropy[n_den_l][n_max], 0.0, beta)
        new_soundspeed= interpolation(soundspeed[n_den_s][n_max], soundspeed[n_den_s][n_max], soundspeed[n_den_l][n_max], soundspeed[n_den_l][n_max], 0.0, beta)        
        #print('a')
    else:
        #print('b')
        temperature_index = temperature_list.index(temperature[temperature>temp_out][0])
        
    
        if temperature_index == 1:
            new_press= interpolation(pressure[n_den_s][0], pressure[n_den_s][0], pressure[n_den_l][0], pressure[n_den_l][0], 0.0, beta)
            new_entropy= interpolation(entropy[n_den_s][0], entropy[n_den_s][0], entropy[n_den_l][0], entropy[n_den_l][0], 0.0, beta)
            new_soundspeed= interpolation(soundspeed[n_den_s][0], soundspeed[n_den_s][0], soundspeed[n_den_l][0], soundspeed[n_den_l][0], 0.0, beta)      
        
        else:

            alpha= (temp_out-temperature[temperature_index-1])/(temperature[temperature_index]-temperature[temperature_index-1])
        
            new_press= interpolation(pressure[n_den_s][temperature_index-1], pressure[n_den_s][temperature_index], pressure[n_den_l][temperature_index-1], pressure[n_den_l][temperature_index], alpha, beta)
            new_entropy= interpolation(entropy[n_den_s][temperature_index-1], entropy[n_den_s][temperature_index], entropy[n_den_l][temperature_index-1], entropy[n_den_l][temperature_index], alpha, beta)
            new_soundspeed= interpolation(soundspeed[n_den_s][temperature_index-1], soundspeed[n_den_s][temperature_index], soundspeed[n_den_l][temperature_index-1], soundspeed[n_den_l][temperature_index], alpha, beta)
    
    
    return temp_out, new_press, new_soundspeed, new_entropy


density_list = list(density)
energy_list = list(energy)
temperature_list = list(temperature)




h = open('test4.txt','w')


# 1D interpolation & extrapolation (linear)
for m in range(0,nr):
    for n in range(0, nu):
        new_temperature[m][n], new_pressure[m][n], new_soundspeed[m][n], new_entropy[m][n] = rho_u(energy_list, density_list, temperature_list, new_energy[n], density[m])

        h.write("%15.8E %15.8E %15.8E  %15.8E  %15.8E  %15.8E  %i %i \n" % (density[m],  new_temperature[m][n], new_energy[n], new_pressure[m][n], new_soundspeed[m][n],  new_entropy[m][n], m,n))
 
#sys.exit()


#new_pressure = new_pressure.clip(min=0.0)







    
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
h.write("%s %i %i %s \n" % ('# ', nr, nu , ':Grid numbers for density and internal energy'))
h.write('# Density (km/m3), Internal energy (kJ/kg), Temperature (K), Pressure (Pa), Sound speed (m/s), Entropy (J/K/kg) \n')


for m in range(0,nr):
    for n in range(0,nu):
        h.write("%15.8E %15.8E %15.8E  %15.8E %15.8E %15.8E \n" % (density[m], new_energy[n],  new_temperature[m][n], new_pressure[m][n], new_soundspeed[m][n], new_entropy[m][n]))
