import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib.gridspec import GridSpec
import gc

round_number=4
timestep1=int(input("Enter initial output number: "))
timestep2=int(input("Enter final output number: "))
ncores=int(input("Enter the number of cores: "))
inputparameter=str(raw_input("Enter Visualization Parameter (Density, Energy, Pressure, Entropy, Temperature, None): "))

if inputparameter!="None":

    center=float(input("Enter the middle % of {} data to visualize: ".format(inputparameter)))
    path=str(raw_input("Enter absolute file path: "))
    outputpath=str(raw_input("Enter absolute snapshot output path: "))
    print

    if inputparameter == "Density":
        parameter=9
        normfactor=1000
        units="1e3 kg/m^3"
    elif inputparameter == "Energy":
        parameter=10
        normfactor=1000000
        units="1E6 J"
    elif inputparameter == "Pressure":
        parameter=11
        normfactor=100000000000
        units="100e9 Pa"
    elif inputparameter == "Entropy":
        parameter=13
        normfactor=1000
        units="1e3 J/K"
    elif inputparameter == "Temperature":
        parameter=14
        normfactor=1000
        units="1e3 K"

    def findmax(timestep1,timestep2,ncores,parameter,center,path):

        valuelist=[]
        for j in range(timestep1,timestep2+1):
            otherfilename="results."+'%05d'%j
            print("Finding min and max {}; reading from {}".format(inputparameter,otherfilename))
            for i in range(ncores):
                filename=path+"results."+'%05d'%j+"_"+'%05d'%ncores+"_"+'%05d'%i+".dat"
                with open (filename,"r") as file:
                    file.readline()
                    file.readline()
                    for line in file:
                        elements=line.split()
                        value=round(float(elements[parameter])/normfactor,round_number)
                        valuelist.append(value)
        lowerper=(100-center)/2
        upperper=center+((100-center)/2)
        lower=np.percentile(valuelist,lowerper)
        upper=np.percentile(valuelist,upperper)
        print("\nThe maximum ({} percentile) normalized value of {} is {}".format(upperper,inputparameter,upper))
        print("The minimum ({} percentile) normalized value of {} is {}\n".format(lowerper,inputparameter,lower))
        del valuelist
        gc.collect()
        return upper,lower

    upper,lower=findmax(timestep1,timestep2,ncores,parameter,center,path)
    for j in range(timestep1,timestep2+1):
        xlist=[]
        ylist=[]
        zlist=[]
        valuelist=[]
        r0=1e6
        N=0

        for i in range(ncores):
            filename=path+"results."+'%05d'%j+"_"+'%05d'%ncores+"_"+'%05d'%i+".dat"
            with open (filename, "r") as file:
                timevalue=str(file.readline().strip())
                num_of_particles=int(file.readline().strip())
                for line in file:
                    elements=line.split()
                    xx=round(float(elements[3])/r0,round_number)
                    yy=round(float(elements[4])/r0,round_number)
                    zz=round(float(elements[5])/r0,round_number)
                    value=round(float(elements[parameter])/normfactor,round_number)
                    xlist.append(xx)
                    ylist.append(yy)
                    zlist.append(zz)
                    valuelist.append(value)
                N=N+num_of_particles
        fig=plt.figure(figsize=(7,6))
        gs=GridSpec(1,2,width_ratios=[20,0.5])
        ax=fig.add_subplot(gs[0],projection='3d')
        norm=plt.Normalize(lower,upper)
        scatter=ax.scatter(xlist,ylist,zlist,c=valuelist,cmap='inferno',norm=norm,s=0.2,edgecolors='none')
        cbarax=fig.add_subplot(gs[1])
        cbar=fig.colorbar(scatter,cax=cbarax)
        pos=cbarax.get_position()
        cbarax.set_position([pos.x0-0.03,pos.y0+0.17,pos.width,pos.height*0.65])
        cbar.set_label("{} ({})".format(inputparameter,units),labelpad=8)
        ax.set_xlabel("x (1e6 m)")
        ax.set_ylabel("y (1e6 m)")
        ax.set_zlabel("z (1e6 m)")
        ticks=np.linspace(-20,20,9)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_zticks(ticks)
        ax.set_xlim(-20,20)
        ax.set_ylim(-20,20)
        ax.set_zlim(-20,20)
        title=ax.set_title("Visualization at t={} s;\nN={} Particles; Middle {}% of {} Data".format(timevalue,N,center,inputparameter))
        titlepos=title.get_position()
        title.set_position((titlepos[0],titlepos[1]+0.05))
        name=outputpath+"snap_"+'%05d'%j+".png"
        plt.savefig(name,dpi=300)
        print("Outputted {}".format(name))
    print

else:

    path=str(raw_input("Enter absolute file path: "))
    outputpath=str(raw_input("Enter snapshot output path: "))
    print

    for j in range(timestep1,timestep2+1):

        xlist=[]
        ylist=[]
        zlist=[]
        r0=1e6
        round_number=4
        N=0

        for i in range(ncores):
            filename=path+"results."+'%05d'%j+"_"+'%05d'%ncores+"_"+'%05d'%i+".dat"
            with open (filename, "r") as file:
                timevalue=str(file.readline().strip())
                num_of_particles=int(file.readline().strip())
                for line in file:
                    elements=line.split()
                    xx = round(float(elements[3])/r0,round_number)
                    yy = round(float(elements[4])/r0,round_number)
                    zz = round(float(elements[5])/r0,round_number)
                    xlist.append(xx)
                    ylist.append(yy)
                    zlist.append(zz)
                N=N+num_of_particles
        fig=plt.figure(figsize=(7,6))
        ax=fig.add_subplot(111,projection='3d')
        ax.scatter(xlist,ylist,zlist,s=0.2,edgecolors='none')
        ax.set_xlabel("x (1e6 m)")
        ax.set_xlabel("y (1e6 m)")
        ax.set_xlabel("z (1e6 m)")
        ticks=np.linspace(-20,20,9)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_zticks(ticks)
        ax.set_xlim(-20,20)
        ax.set_ylim(-20,20)
        ax.set_zlim(-20,20)
        title=ax.set_title("Visualization at t={} s;\nN={} Particles".format(timevalue,N))
        titlepos=title.get_position()
        title.set_position((titlepos[0],titlepos[1]+0.05))
        name=outputpath+"snap_"+'%05d'%j+".png"
        plt.savefig(name,dpi=300)
        print("Outputted {}".format(name))
    print
