import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

#USER INPUTS ---------------------------------------------

path = 
outputpath =
ncores = 
outputnumber1 = 
outputnumber2 =
axesdim = 
axesscale = 
axes = 
elevation = 
azimuth = 
background = 
particlesize = 

parameter = 
minimum = 
maximum = 
colormap = 

tarmantlecolor = 
tarcorecolor = 
impmantlecolor = 
impcorecolor = 

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

for j in range(outputnumber1,outputnumber2+1):

    xlist=[]
    ylist=[]
    zlist=[]

    xtarmant=[]
    ytarmant=[]
    ztarmant=[]
    tarmantvalue=[]

    xtarcore=[]
    ytarcore=[]
    ztarcore=[]
    tarcorevalue=[]

    ximpmant=[]
    yimpmant=[]
    zimpmant=[]
    impmantvalue=[]

    ximpcore=[]
    yimpcore=[]
    zimpcore=[]
    impcorevalue=[]


    for i in range(ncores):

        filename=f"{path}/results.{j:05d}_{ncores:05d}_{i:05d}.dat"

        with open (filename, "r") as file:
            timevalue=float(file.readline().strip())
            file.readline()
        
            for line in file:
                elements=line.split()
                tag=int(elements[1])
                xx = float(elements[3])/axesscale
                yy = float(elements[4])/axesscale
                zz = float(elements[5])/axesscale
                value = float(elements[index])/normfactor

                if tag == 0:
                    xtarmant.append(xx)
                    ytarmant.append(yy)
                    ztarmant.append(zz)
                    tarmantvalue.append(value)

                if tag == 1:
                    xtarcore.append(xx)
                    ytarcore.append(yy)
                    ztarcore.append(zz)
                    tarcorevalue.append(value)

                if tag == 2:

                    ximpmant.append(xx)
                    yimpmant.append(yy)
                    zimpmant.append(zz)
                    impmantvalue.append(value)

                if tag == 3:

                    ximpcore.append(xx)
                    yimpcore.append(yy)
                    zimpcore.append(zz)
                    impcorevalue.append(value)

    if background == "Black":
        plt.style.use('dark_background')

    fig=plt.figure(figsize=(7,6))
    gs=GridSpec(1,2,width_ratios=[20,0.5],figure=fig)
    ax=fig.add_subplot(gs[0],projection='3d')
    ax.view_init(elev=elevation, azim=azimuth, roll=0)
    norm=plt.Normalize(minimum/normfactor, maximum/normfactor)
    
    if tarmantlecolor == "cmap":
        scatter=ax.scatter(xtarmant,ytarmant,ztarmant,c=tarmantvalue,norm=norm,cmap=colormap,marker='.',s=particlesize)
    else:
        ax.scatter(xtarmant,ytarmant,ztarmant,c=tarmantlecolor,marker='.',s=particlesize)

    if tarcorecolor == "cmap":
        scatter=ax.scatter(xtarcore,ytarcore,ztarcore,c=tarcorevalue,cmap=colormap,norm=norm,marker='.',s=particlesize)
    else:
        ax.scatter(xtarcore,ytarcore,ztarcore,c=tarcorecolor,marker='.',s=particlesize)

    if impmantlecolor == "cmap":
        scatter=ax.scatter(ximpmant,yimpmant,zimpmant,c=impmantvalue,cmap=colormap,norm=norm,marker='.',s=particlesize)
    else:
        ax.scatter(ximpmant,yimpmant,zimpmant,c=impmantlecolor,marker='.',s=particlesize)
    
    if impcorecolor == "cmap":
        scatter=ax.scatter(ximpcore,yimpcore,zimpcore,c=impcorevalue,cmap=colormap,norm=norm,marker='.',s=particlesize)
    else:
        ax.scatter(ximpcore,yimpcore,zimpcore,c=impcorecolor,marker='.',s=particlesize)

    if tarmantlecolor == 'cmap' or tarcorecolor == 'cmap' or impmantlecolor == 'cmap' or impcorecolor == 'cmap':

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
            cbarax.set_position([pos.x0,pos.y0+0.13,pos.width,pos.height*0.65])
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

    if tarmantlecolor == 'cmap' or tarcorecolor == 'cmap' or impmantlecolor == 'cmap' or impcorecolor == 'cmap':
        title = ax.set_title(f"       Visualization of {parameter} at t={timevalue:.2e}")
    else:
        title = ax.set_title(f"Visualization at t={timevalue:.2e}")

    pos = ax.get_position()
    title.set_position([0.5, 1.02])
    ax.set_box_aspect([1,1,1])

    plt.savefig(f"{outputpath}/snap_{j:05d}.png", dpi=300)
    print(f"Outputted snap_{j:05d}.png")
    #plt.show()
    plt.close()



