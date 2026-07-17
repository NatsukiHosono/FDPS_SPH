import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

#USER INPUTS_________________________________________________________________

path = '/scratch/mnakaji2/FDPS_SPH_CMAP/FDPS_SPH/output'
outputpath = '/home/mnakaji2'
outputnumber1 = 0
outputnumber2 = 5
ncores = 80

centering = True
axesscale = 1e6
axesdim = 20
axes = False
background = 'Black'
thickness = 1
particlesize = 1


# -- mantle property
parameter = 'Temperature'
minimum = 1000
maximum = 4000
colormap = 'afmhot'

# -- core property

parameter2 = 'Temperature'
minimum2 = 2000
maximum2 = 4000
colormap2 = 'PiYG'

tarmantlecolor = 'cmap'
tarcorecolor = 'cmap'
impmantlecolor = 'cmap'
impcorecolor = 'cmap'

#____________________________________________________________________________

def normalize(parameter):    # Determines parameter normalization factor and units

    if parameter=='Density':

        index=9
        normfactor=1000
        units='1e3 $\mathrm{kg/m^3}$'

    elif parameter=='Energy':

        index=10
        normfactor=1000000
        units='1E6 $\mathrm{J}$'

    elif parameter=='Pressure':

        index=11
        normfactor=1000000000
        units='$\mathrm{GPa}$'

    elif parameter=='Entropy':

        index=13
        normfactor=1000
        units='1e3 $\mathrm{J/K}$'

    elif parameter=='Temperature':

        index=14
        normfactor=1000
        units='1e3 $\mathrm{K}$'

    return index,normfactor,units



def read_data(index,normfactor,j):    # Reads and stores particle data

    zmin=-thickness/2
    zmax=thickness/2

    xtarmant=[]
    ytarmant=[]
    tarmantvalue=[]

    xtarcore=[]
    ytarcore=[]
    tarcorevalue=[]

    ximpmant=[]
    yimpmant=[]
    impmantvalue=[]

    ximpcore=[]
    yimpcore=[]
    impcorevalue=[]

    masslist=[]

    for i in range(ncores):

        filename=f'{path}/results.{j:05d}_{ncores:05d}_{i:05d}.dat'

        with open (filename, 'r') as file:

            timevalue=float(file.readline().strip())
            file.readline()
        
            for line in file:

                elements=line.split()
                tag=int(elements[1])
                m=float(elements[2])
                xx=float(elements[3])/axesscale
                yy=float(elements[4])/axesscale
                zz=float(elements[5])/axesscale
                value=float(elements[index])/normfactor

                if zmin < zz < zmax:

                    if tag == 0:

                        xtarmant.append(xx)
                        ytarmant.append(yy)
                        tarmantvalue.append(value)

                    if tag == 1:

                        xtarcore.append(xx)
                        ytarcore.append(yy)
                        tarcorevalue.append(value)
                        masslist.append(m)

                    if tag == 2:

                        ximpmant.append(xx)
                        yimpmant.append(yy)
                        impmantvalue.append(value)

                    if tag == 3:

                        ximpcore.append(xx)
                        yimpcore.append(yy)
                        impcorevalue.append(value)

    lists=[np.array(xtarmant),np.array(ytarmant),np.array(tarmantvalue),
           np.array(xtarcore),np.array(ytarcore),np.array(tarcorevalue),
           np.array(ximpmant),np.array(yimpmant),np.array(impmantvalue),
           np.array(ximpcore),np.array(yimpcore),np.array(impcorevalue)]
    
    return lists,timevalue,masslist



def center_particles(lists,masslist):    # Transforms to center of mass frame

    masslist=np.array(masslist)
    
    xtarmant,ytarmant,tarmantvalue=lists[0],lists[1],lists[2]
    xtarcore,ytarcore,tarcorevalue=lists[3],lists[4],lists[5]
    ximpmant,yimpmant,impmantvalue=lists[6],lists[7],lists[8]
    ximpcore,yimpcore,impcorevalue=lists[9],lists[10],lists[11]

    xcm=np.sum(xtarcore*masslist)/np.sum(masslist)
    ycm=np.sum(ytarcore*masslist)/np.sum(masslist)

    xtarmant,ytarmant=xtarmant-xcm,ytarmant-ycm
    xtarcore,ytarcore=xtarcore-xcm,ytarcore-ycm
    ximpmant,yimpmant=ximpmant-xcm,yimpmant-ycm
    ximpcore,yimpcore=ximpcore-xcm,yimpcore-ycm

    lists=lists.clear()

    lists=[xtarmant,ytarmant,tarmantvalue,
           xtarcore,ytarcore,tarcorevalue,
           ximpmant,yimpmant,impmantvalue,
           ximpcore,yimpcore,impcorevalue]
    
    return lists



def plot(lists,normfactor,units,timevalue,j):    # Plots the cross section

    xtarmant,ytarmant,tarmantvalue=lists[0],lists[1],lists[2]
    xtarcore,ytarcore,tarcorevalue=lists[3],lists[4],lists[5]
    ximpmant,yimpmant,impmantvalue=lists[6],lists[7],lists[8]
    ximpcore,yimpcore,impcorevalue=lists[9],lists[10],lists[11]

    if background == 'Black':

        plt.style.use('dark_background')

    fig,ax=plt.subplots(figsize=(2112/300,1808/300))
    gs=GridSpec(1,2,width_ratios=[20,0.5],figure=fig)
    norm=plt.Normalize(minimum/normfactor, maximum/normfactor)
    norm2=plt.Normalize(minimum2/normfactor, maximum2/normfactor)
    
    if tarmantlecolor == 'cmap':
        scatter=ax.scatter(xtarmant,ytarmant,c=tarmantvalue,norm=norm,cmap=colormap,marker='.',s=particlesize)
    else:
        ax.scatter(xtarmant,ytarmant,c=tarmantlecolor,marker='.',s=particlesize)

    if tarcorecolor == 'cmap':
        scatter2=ax.scatter(xtarcore,ytarcore,c=tarcorevalue,cmap=colormap2,norm=norm2,marker='.',s=particlesize)
    else:
        ax.scatter(xtarcore,ytarcore,c=tarcorecolor,marker='.',s=particlesize)

    if impmantlecolor == 'cmap':
        scatter=ax.scatter(ximpmant,yimpmant,c=impmantvalue,cmap=colormap,norm=norm,marker='.',s=particlesize)
    else:
        ax.scatter(ximpmant,yimpmant,c=impmantlecolor,marker='.',s=particlesize)
    
    if impcorecolor == 'cmap':
        scatter2=ax.scatter(ximpcore,yimpcore,c=impcorevalue,cmap=colormap2,norm=norm2,marker='.',s=particlesize)
    else:
        ax.scatter(ximpcore,yimpcore,c=impcorecolor,marker='.',s=particlesize)

    if tarmantlecolor == 'cmap' or tarcorecolor == 'cmap' or impmantlecolor == 'cmap' or impcorecolor == 'cmap':

        if axes == False:
                
            cbar_width = 0.3
            cbar_height = 0.025
            cbar_bottom = 0.15
            
            cbarax1 = fig.add_axes([0.17, cbar_bottom, cbar_width, cbar_height])
            cbar1 = fig.colorbar(scatter, cax=cbarax1, orientation='horizontal')
            cbar1.set_label(f"Mantle {parameter} ({units})", labelpad=8)

            cbarax2 = fig.add_axes([0.55, cbar_bottom, cbar_width, cbar_height])
            cbar2 = fig.colorbar(scatter2, cax=cbarax2, orientation='horizontal')
            cbar2.set_label(f"Core {parameter2} ({units})", labelpad=8)
            

        else:

            cbarax1 = fig.add_axes([0.84, 0.25, 0.025, 0.5])
            cbar1 = fig.colorbar(scatter, cax=cbarax1)
            cbar1.set_label(f"Mantle {parameter} ({units})", labelpad=8)

            cbarax2 = fig.add_axes([0.92, 0.25, 0.025, 0.5])
            cbar2 = fig.colorbar(scatter2, cax=cbarax2)
            cbar2.set_label(f"Core {parameter2} ({units})", labelpad=8)

            

    ax.set_xlabel(f'x ({axesscale:.0e} m)')
    ax.set_ylabel(f'y ({axesscale:.0e} m)')

    ticks=np.linspace(-axesdim/2,axesdim/2,5)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    ax.set_xlim(-axesdim/2,axesdim/2)
    ax.set_ylim(-axesdim/2,axesdim/2)

    if axes == False:

        ax.set_axis_off()

    if tarmantlecolor == 'cmap' or tarcorecolor == 'cmap' or impmantlecolor == 'cmap' or impcorecolor == 'cmap':

        title = ax.set_title(f'Cross Section of {parameter} at t={timevalue/3600:0.2f} h')

    else:

        title = ax.set_title(f'Cross Section at t={timevalue/3600:0.2f} h')

    pos = ax.get_position()
    title.set_position([0.5, 1.02])
    ax.set_aspect('equal')

    plt.savefig(f'{outputpath}/CS_{j:05d}.png', dpi=300)

    print(f'    Outputted CS_{j:05d}.png')

    plt.close()



def main():
    
    print('\nGenerating cross sections...\n')

    for j in range(outputnumber1,outputnumber2+1):

        index,normfactor,units=normalize(parameter)
        lists,timevalue,masslist=read_data(index,normfactor,j)

        if centering == True:

            lists=center_particles(lists,masslist)

        plot(lists,normfactor,units,timevalue,j)

    print()



main()
