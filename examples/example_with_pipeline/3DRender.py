import numpy as np
import pyvista as pv
import vtk
vtk.vtkObject.GlobalWarningDisplayOff()

#USER INPUTS_________________________________________________________________

path = r"C:\Users\nakaj\OneDrive\Desktop\2025Work\MeltScalingCollisions\FinalCollisions\0.3\60Deg\\"
outputpath = path
ncores = 150
outputnumber1 = 740
outputnumber2 = 740

centering = True
axesscale = 1e6
cameraposition = (110,110,0)
axes = False
azimuth = 0
elevation = 20
background = 'Black'
particlesize = 3

parameter = 'Energy'
minimum = 0
maximum = 30e6
colormap = 'magma'

tarmantlecolor = 'cmap'
tarcorecolor = 'cmap'
impmantlecolor = 'cmap'
impcorecolor = 'cmap'

#____________________________________________________________________________

def normalize(parameter):    # Determines parameter normalization factor and units

    if parameter=='Density':

        index=9
        normfactor=1000
        units='1e3 kg/m³$'

    elif parameter=='Energy':

        index=10
        normfactor=1000000
        units='1E6 J'

    elif parameter=='Pressure':

        index=11
        normfactor=1000000000
        units='GPa'

    elif parameter=='Entropy':

        index=13
        normfactor=1000
        units='1e3 J/K'

    elif parameter=='Temperature':

        index=14
        normfactor=1000
        units='1e3 K'

    return index,normfactor,units



def read_data(index,normfactor,j):    # Reads and stores particle data

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
                    masslist.append(m)

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

    lists=[np.array(xtarmant),np.array(ytarmant),np.array(ztarmant),np.array(tarmantvalue),
        np.array(xtarcore),np.array(ytarcore),np.array(ztarcore),np.array(tarcorevalue),
        np.array(ximpmant),np.array(yimpmant),np.array(zimpmant),np.array(impmantvalue),
        np.array(ximpcore),np.array(yimpcore),np.array(zimpcore),np.array(impcorevalue)]
    
    return lists,timevalue,masslist         



def center_particles(lists,masslist):    # Transforms to center of mass frame

    masslist=np.array(masslist)
    
    xtarmant,ytarmant,ztarmant,tarmantvalue=lists[0],lists[1],lists[2],lists[3]
    xtarcore,ytarcore,ztarcore,tarcorevalue=lists[4],lists[5],lists[6],lists[7]
    ximpmant,yimpmant,zimpmant,impmantvalue=lists[8],lists[9],lists[10],lists[11]
    ximpcore,yimpcore,zimpcore,impcorevalue=lists[12],lists[13],lists[14],lists[15]

    xcm=np.sum(xtarcore*masslist)/np.sum(masslist)
    ycm=np.sum(ytarcore*masslist)/np.sum(masslist)
    zcm=np.sum(ztarcore*masslist)/np.sum(masslist)

    xtarmant,ytarmant,ztarmant=xtarmant-xcm,ytarmant-ycm,ztarmant-zcm
    xtarcore,ytarcore,ztarcore=xtarcore-xcm,ytarcore-ycm,ztarcore-zcm
    ximpmant,yimpmant,zimpmant=ximpmant-xcm,yimpmant-ycm,zimpmant-zcm
    ximpcore,yimpcore,zimpcore=ximpcore-xcm,yimpcore-ycm,zimpcore-zcm

    lists=[xtarmant,ytarmant,ztarmant,tarmantvalue,
           xtarcore,ytarcore,ztarcore,tarcorevalue,
           ximpmant,yimpmant,zimpmant,impmantvalue,
           ximpcore,yimpcore,zimpcore,impcorevalue]
    
    return lists



def plot(lists,normfactor,units,timevalue,j):    # Plots the 3D Render

    xtarmant, ytarmant, ztarmant, tarmantvalue = lists[0], lists[1], lists[2], lists[3]
    xtarcore, ytarcore, ztarcore, tarcorevalue = lists[4], lists[5], lists[6], lists[7]
    ximpmant, yimpmant, zimpmant, impmantvalue = lists[8], lists[9], lists[10], lists[11]
    ximpcore, yimpcore, zimpcore, impcorevalue = lists[12], lists[13], lists[14], lists[15]

    p = pv.Plotter(off_screen=True, window_size=(2112, 1808)) 
    pv.global_theme.font.family = 'courier'

    cmap_min = minimum / normfactor
    cmap_max = maximum / normfactor

    width=0.4
    height=0.06
    position_x=0.5-width/2

    if background == 'Black':

        textcolor = 'white'

    if background == 'White':

        textcolor = 'black'

    if tarmantlecolor == 'cmap':

        points = np.column_stack([xtarmant, ytarmant, ztarmant])
        cloud = pv.PolyData(points)
        cloud['values'] = tarmantvalue
        p.add_points(cloud,scalars='values',cmap=colormap,clim=[cmap_min, cmap_max],
                     point_size=particlesize,render_points_as_spheres=True,lighting=False,show_scalar_bar=False)
        p.add_scalar_bar(title=f"{parameter} ({units})\n",n_labels=6,vertical=False,color=textcolor,
                         title_font_size=30,label_font_size=26,position_x=position_x,width=width,height=height)

    else:

        points = np.column_stack([xtarmant, ytarmant, ztarmant])
        p.add_points(points, color=tarmantlecolor, point_size=particlesize,render_points_as_spheres=True,lighting=False,show_scalar_bar=False)


    if tarcorecolor == 'cmap':

        points = np.column_stack([xtarcore, ytarcore, ztarcore])
        cloud = pv.PolyData(points)
        cloud['values'] = tarcorevalue
        p.add_points(cloud, scalars='values', cmap=colormap, clim=[cmap_min, cmap_max],
                     point_size=particlesize, render_points_as_spheres=True,lighting=False,show_scalar_bar=False)
        p.remove_scalar_bar()
        p.add_scalar_bar(title=f"{parameter} ({units})\n",n_labels=6,vertical=False,color=textcolor,
                         title_font_size=30,label_font_size=26,position_x=position_x,width=width,height=height)
        
    else:

        points = np.column_stack([xtarcore, ytarcore, ztarcore])
        p.add_points(points, color=tarcorecolor, point_size=particlesize,render_points_as_spheres=True,lighting=False,show_scalar_bar=False)

   
    if impmantlecolor == 'cmap':

        points = np.column_stack([ximpmant, yimpmant, zimpmant])
        cloud = pv.PolyData(points)
        cloud['values'] = impmantvalue
        p.add_points(cloud, scalars='values', cmap=colormap, clim=[cmap_min, cmap_max],
                     point_size=particlesize,render_points_as_spheres=True,lighting=False,show_scalar_bar=False)
        p.remove_scalar_bar()
        p.add_scalar_bar(title=f"{parameter} ({units})\n",n_labels=6,vertical=False,color=textcolor,
                         title_font_size=30,label_font_size=26,position_x=position_x,width=width,height=height)
        
    else:

        points = np.column_stack([ximpmant, yimpmant, zimpmant])
        p.add_points(points, color=impmantlecolor, point_size=particlesize,render_points_as_spheres=True,lighting=False,show_scalar_bar=False)


    if impcorecolor == 'cmap':

        points = np.column_stack([ximpcore, yimpcore, zimpcore])
        cloud = pv.PolyData(points)
        cloud['values'] = impcorevalue
        p.add_points(cloud, scalars='values', cmap=colormap, clim=[cmap_min, cmap_max],
                     point_size=particlesize,render_points_as_spheres=True,lighting=False,show_scalar_bar=False)
        p.remove_scalar_bar()
        p.add_scalar_bar(title=f"{parameter} ({units})\n",n_labels=6,vertical=False,color=textcolor,
                         title_font_size=30,label_font_size=26,position_x=position_x,width=width,height=height)
        
    else:

        points = np.column_stack([ximpcore, yimpcore, zimpcore])
        p.add_points(points, color=impcorecolor,point_size=particlesize,render_points_as_spheres=True,lighting=False,show_scalar_bar=False)

    if axes:

        p.add_axes(color=textcolor,line_width=2,labels_off=False)

    if background == 'Black':

        p.set_background('black')
        
    p.add_title(title=f'3D Render of {parameter} at {timevalue/3600:0.2f} h',font_size= 18,color=textcolor)
    p.camera_position = [cameraposition,(0,0,0),(0,0,1)]  
    p.camera.Azimuth(azimuth)
    p.camera.Elevation(elevation)

    p.show(screenshot=f"{outputpath}/3D_{j:05d}.png")

    p.close()

    print(f'    Outputted 3D_{j:05d}.png')



def main():
    
    print('\nGenerating 3D renders...\n')

    for j in range(outputnumber1,outputnumber2+1):

        index,normfactor,units=normalize(parameter)
        lists,timevalue,masslist=read_data(index,normfactor,j)

        if centering == True:

            lists=center_particles(lists,masslist)

        plot(lists,normfactor,units,timevalue,j)

    print()



main()