import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

#USER INPUTS ---------------------------------------------

path =
outputpath =
ncores = 
timestep1 =
timestep2 =
parameter =
minimum =
maximum =
axesdim =
axesscale =
axes =
background =

#---------------------------------------------------------

if parameter=="Density":
    index=9
    normfactor=1000
    units="1e3 $\mathrm{kg/m^3}$"
elif parameter=="Energy":
    index=10
    normfactor=1000000
    units="1E6 $\mathrm{J}$"
elif parameter=="Pressure":
    index=11
    normfactor=100000000000
    units="100e9 $\mathrm{Pa}$"
elif parameter=="Entropy":
    index=13
    normfactor=1000
    units="1e3 $\mathrm{J/K}$"
elif parameter=="Temperature":
    index=14
    normfactor=1000
    units="1e3 $\mathrm{K}$"

if parameter != "None":

    for i in range(timestep1,timestep2+1):

        xlist=[]
        ylist=[]
        zlist=[]
        valuelist=[]

        for j in range(ncores):

            filename=fr"{path}/results.{i:05d}_{ncores:05d}_{j:05d}.dat"

            with open (filename, "r") as file:
                timevalue=float(file.readline().strip())
                file.readline()
            
                for line in file:
                    elements=line.split()
                    x=float(elements[3])/axesscale
                    y=float(elements[4])/axesscale
                    z=float(elements[5])/axesscale
                    value=float(elements[index])/normfactor
                    valuelist.append(value)
                    xlist.append(x)
                    ylist.append(y)
                    zlist.append(z)

        if background == "Black":
            plt.style.use('dark_background')

        fig=plt.figure(figsize=(7,6))
        gs=GridSpec(1,2,width_ratios=[20,0.5],figure=fig)
        ax=fig.add_subplot(gs[0],projection='3d')

        norm=plt.Normalize(minimum/normfactor, maximum/normfactor)
        scatter=ax.scatter(xlist,ylist,zlist,c=valuelist,cmap='inferno',norm=norm,marker='.',s=0.07)
        
        if axes == "No":
            cbarax=fig.add_subplot(gs[1])
            cbar=fig.colorbar(scatter,cax=cbarax,orientation='horizontal')
            cbar_width = 0.4
            cbar_height = 0.025
            cbar_bottom = 0.15
            cbar_left = (1 - cbar_width) / 2
            cbarax.set_position([cbar_left, cbar_bottom, cbar_width, cbar_height])
            cbar.set_label(f"{parameter} ({units})", labelpad=8)

        else:

            cbarax=fig.add_subplot(gs[1])
            cbar=fig.colorbar(scatter,cax=cbarax)
            pos=cbarax.get_position()
            cbarax.set_position([pos.x0,pos.y0+0.17,pos.width,pos.height*0.65])
            cbar.set_label(f"{parameter} ({units})", labelpad=8)

        ax.set_xlabel(f"x ({axesscale:.0e} m)")
        ax.set_ylabel(f"y ({axesscale:.0e} m)")
        ax.set_zlabel(f"z ({axesscale:.0e} m)")

        ticks=np.linspace(-axesdim/2,axesdim/2,5)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_zticks(ticks)

        ax.set_xlim(-axesdim/2,axesdim/2)
        ax.set_ylim(-axesdim/2,axesdim/2)
        ax.set_zlim(-axesdim/2,axesdim/2)

        if axes == "No":
            ax.set_axis_off()

        title = ax.set_title(f"       Visualization of {parameter} at t={timevalue:.2e}")
        pos = ax.get_position()
        title.set_position([0.5, 1.02])
        ax.set_box_aspect([1,1,1])
        print(f"Outputted snap_{i:05d}.png")
        plt.savefig(f"{outputpath}/snap_{i:05d}.png", dpi=300)
        plt.close()

else:

    for j in range(timestep1,timestep2+1):

        xlist=[]
        ylist=[]
        zlist=[]

        for i in range(ncores):

            filename=f"{path}/results.{j:05d}_{ncores:05d}_{i:05d}.dat"

            with open (filename, "r") as file:
                timevalue=float(file.readline().strip())
                file.readline()
            
                for line in file:
                    elements=line.split()
                    xx = float(elements[3])/axesscale
                    yy = float(elements[4])/axesscale
                    zz = float(elements[5])/axesscale

                    xlist.append(xx)
                    ylist.append(yy)
                    zlist.append(zz)

        if background == "Black":
            plt.style.use('dark_background')

        fig=plt.figure(figsize=(5,5))
        ax=fig.add_subplot(111,projection='3d')
        ax.scatter(xlist,ylist,zlist,marker='.',s=0.07)  

        ax.set_xlabel(f"x ({axesscale:.0e} m)")
        ax.set_ylabel(f"y ({axesscale:.0e} m)")
        ax.set_zlabel(f"z ({axesscale:.0e} m)")

        ticks=np.linspace(-axesdim/2,axesdim/2,5)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_zticks(ticks)

        ax.set_xlim(-axesdim/2,axesdim/2)
        ax.set_ylim(-axesdim/2,axesdim/2)
        ax.set_zlim(-axesdim/2,axesdim/2)

        if axes == "No":
            ax.set_axis_off()

        ax.set_title(f"Visualization at t={timevalue:.2e}")

        ax.set_box_aspect([1,1,1])

        print(f"Outputted snap_{j:05d}.png")
        plt.savefig(f"{outputpath}/snap_{j:05d}.png", dpi=300)
        plt.close()



