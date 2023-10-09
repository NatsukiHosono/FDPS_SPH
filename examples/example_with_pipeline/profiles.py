import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib, sys

if len(sys.argv) != 4:
    print "Usage inputs: <path> <time [s]> <ncores>"
    exit(1)

path = sys.argv[1]
timestep = np.int(np.float(sys.argv[2]))/100
n = np.int(np.float(sys.argv[3]))

def readfile(path,step,ncores):
    frames=[pd.read_csv(path+"results."+'%05d'%step+"_"+'%05d'%ncores+"_"'%05d'%f+".dat",skiprows=2,header=None,sep='\t') for f in range(ncores)]
    b=pd.concat(frames,ignore_index=True).sort_values(by=[0])
    b.rename(columns={0:'id', 1:'tag',  2:'mass', 3:'x', 4:'y', 5:'z', 6:'x vel', 7:'y vel', 8:'z vel', 9:'dens', 10:'int ener', 11:'pres', 12:'pot ener', 13:'entr', 14:'temp'}, inplace=True)
    b['r']=(b['x']**2+b['y']**2+b['z']**2)**0.5/1000
    b['v']=(b['x vel']**2+b['y vel']**2+b['z vel']**2)**0.5
    return b

def profileplot(b,step):
    fig, axs = plt.subplots(2,2)
    fig.set_size_inches(8.6,8.6)
    fig.tight_layout()
    for i in range (2):
        for j in range(2):
	     axs[i,j].grid(True)
             axs[i,j].set_xlabel('r [km]')
    axs[0,0].plot(b['r'],b['dens']/1000,'.',color='#a6bddb')
    axs[0,0].set_ylabel('density [g/cm3]')
    axs[0,1].plot(b['r'],b['temp'],'.',color='#a6bddb')
    axs[0,1].set_ylabel('temperature [K]')
    axs[1,0].plot(b['r'],b['pres']/1000000000,'.',color='#a6bddb')
    axs[1,0].set_ylabel('pressure [GPa]')
    axs[1,1].plot(b['r'],b['v'],'.',color='#a6bddb')
    axs[1,1].set_ylabel('velocity [m/s]')
    plt.savefig('profile_'+'%05d'%step+'.png',dpi=300,bbox_inches='tight',facecolor='white',transparent=False)

profileplot(readfile(path,timestep,n),timestep)
