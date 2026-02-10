import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

# USER INPUTS_________________________________________________________________

path = '/scratch/mnakaji2/FDPS_SPH3/FDPS_SPH/output' 
outputpath = '/home/mnakaji2'
outputnumber1 = 0
outputnumber2 = 100
ncores = 40

centering = True
axesscale = 1e6
axesdim = 40
axes = False
background = 'Black'
thickness = 10
particlesize = 1

parameter = 'Temperature'
minimum = 2000
maximum = 8000
colormap = 'afmhot'

tarmantlecolor = 'cmap'
tarcorecolor = 'gray'
impmantlecolor = 'cmap'
impcorecolor = 'gray'

#____________________________________________________________________________

def normalize(parameter):    # Determines parameter normalization factor and units
    if parameter == 'Density':
        index = 9
        normfactor = 1000
        units = '1e3 $\\mathrm{kg/m^3}$'

    elif parameter == 'Energy':
        index = 10
        normfactor = 1000000
        units = '1E6 $\\mathrm{J}$'

    elif parameter == 'Pressure':
        index = 11
        normfactor = 1000000000
        units = '$\\mathrm{GPa}$'

    elif parameter == 'Entropy':
        index = 13
        normfactor = 1000.0
        units = '1e3 $\\mathrm{J/K}$'

    elif parameter == 'Temperature':
        index = 14
        normfactor = 1000
        units = '1e3 $\\mathrm{K}$'

    return index, normfactor, units


def read_data(index, normfactor, j):    # Reads and stores particle data
    zmin = -thickness / 2.0
    zmax = thickness / 2.0

    xtarmant, ytarmant, tarmantvalue = [], [], []
    xtarcore, ytarcore, tarcorevalue = [], [], []
    ximpmant, yimpmant, impmantvalue = [], [], []
    ximpcore, yimpcore, impcorevalue = [], [], []

    masslist = []

    for i in range(ncores):
        filename = '{0}/results.{1:05d}_{2:05d}_{3:05d}.dat'.format(path, j, ncores, i)

        with open(filename, 'r') as file:
            timevalue = float(file.readline().strip())
            file.readline()

            for line in file:
                elements = line.split()
                tag = int(elements[1])
                m = float(elements[2])
                xx = float(elements[3]) / axesscale
                yy = float(elements[4]) / axesscale
                zz = float(elements[5]) / axesscale
                value = float(elements[index]) / normfactor

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

    lists = [np.array(xtarmant), np.array(ytarmant), np.array(tarmantvalue),
             np.array(xtarcore), np.array(ytarcore), np.array(tarcorevalue),
             np.array(ximpmant), np.array(yimpmant), np.array(impmantvalue),
             np.array(ximpcore), np.array(yimpcore), np.array(impcorevalue)]

    return lists, timevalue, masslist


def center_particles(lists, masslist):    # Transforms to center of mass frame
    masslist = np.array(masslist)

    xtarmant, ytarmant, tarmantvalue = lists[0], lists[1], lists[2]
    xtarcore, ytarcore, tarcorevalue = lists[3], lists[4], lists[5]
    ximpmant, yimpmant, impmantvalue = lists[6], lists[7], lists[8]
    ximpcore, yimpcore, impcorevalue = lists[9], lists[10], lists[11]

    xcm = np.sum(xtarcore * masslist) / np.sum(masslist)
    ycm = np.sum(ytarcore * masslist) / np.sum(masslist)

    xtarmant, ytarmant = xtarmant - xcm, ytarmant - ycm
    xtarcore, ytarcore = xtarcore - xcm, ytarcore - ycm
    ximpmant, yimpmant = ximpmant - xcm, yimpmant - ycm
    ximpcore, yimpcore = ximpcore - xcm, yimpcore - ycm

    # IMPORTANT: do NOT do lists = lists.clear() (that returns None)
    lists = [xtarmant, ytarmant, tarmantvalue,
             xtarcore, ytarcore, tarcorevalue,
             ximpmant, yimpmant, impmantvalue,
             ximpcore, yimpcore, impcorevalue]

    return lists


def plot(lists, normfactor, units, timevalue, j):    # Plots the cross section
    xtarmant, ytarmant, tarmantvalue = lists[0], lists[1], lists[2]
    xtarcore, ytarcore, tarcorevalue = lists[3], lists[4], lists[5]
    ximpmant, yimpmant, impmantvalue = lists[6], lists[7], lists[8]
    ximpcore, yimpcore, impcorevalue = lists[9], lists[10], lists[11]

    if background == 'Black':
        plt.style.use('dark_background')

    fig, ax = plt.subplots(figsize=(2112/300.0, 1808/300.0))
    gs = GridSpec(1, 2, width_ratios=[20, 0.5])
    norm = plt.Normalize(minimum / float(normfactor), maximum / float(normfactor))

    scatter = None

    if tarmantlecolor == 'cmap':
        scatter = ax.scatter(xtarmant, ytarmant, c=tarmantvalue, norm=norm,
                             cmap=colormap, marker='.', s=particlesize,edgecolors='none')
    else:
        ax.scatter(xtarmant, ytarmant, c=tarmantlecolor, marker='.', s=particlesize,edgecolors='none')

    if tarcorecolor == 'cmap':
        scatter = ax.scatter(xtarcore, ytarcore, c=tarcorevalue, cmap=colormap,
                             norm=norm, marker='.', s=particlesize, edgecolors='none')
    else:
        ax.scatter(xtarcore, ytarcore, c=tarcorecolor, marker='.', s=particlesize,edgecolors='none')

    if impmantlecolor == 'cmap':
        scatter = ax.scatter(ximpmant, yimpmant, c=impmantvalue, cmap=colormap,
                             norm=norm, marker='.', s=particlesize,edgecolors='none')
    else:
        ax.scatter(ximpmant, yimpmant, c=impmantlecolor, marker='.', s=particlesize, edgecolors='none')

    if impcorecolor == 'cmap':
        scatter = ax.scatter(ximpcore, yimpcore, c=impcorevalue, cmap=colormap,
                             norm=norm, marker='.', s=particlesize,edgecolors='none')
    else:
        ax.scatter(ximpcore, yimpcore, c=impcorecolor, marker='.', s=particlesize,edgecolors='none')

    using_cmap = (tarmantlecolor == 'cmap' or tarcorecolor == 'cmap' or
                  impmantlecolor == 'cmap' or impcorecolor == 'cmap')

    if using_cmap and scatter is not None:
        cbarax = fig.add_subplot(gs[1])

        if not axes:
            cbar = fig.colorbar(scatter, cax=cbarax, orientation='horizontal')
            cbar_width = 0.4
            cbar_height = 0.025
            cbar_bottom = 0.15
            cbar_left = (1 - cbar_width) / 2.0
            cbarax.set_position([cbar_left, cbar_bottom, cbar_width, cbar_height])
            cbar.set_label("{0} ({1})".format(parameter, units), labelpad=8)
        else:
            cbar = fig.colorbar(scatter, cax=cbarax)
            pos = cbarax.get_position()
            cbarax.set_position([pos.x0, pos.y0 + 0.13, pos.width, pos.height * 0.65])
            cbar.set_label("{0} ({1})".format(parameter, units), labelpad=8)

    ax.set_xlabel('x ({0:.0e} m)'.format(axesscale))
    ax.set_ylabel('y ({0:.0e} m)'.format(axesscale))

    ticks = np.linspace(-axesdim/2.0, axesdim/2.0, 5)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    ax.set_xlim(-axesdim/2.0, axesdim/2.0)
    ax.set_ylim(-axesdim/2.0, axesdim/2.0)

    if not axes:
        ax.set_axis_off()

    if using_cmap:
        title = ax.set_title('Cross Section of {0} at t={1:0.2f} h'.format(parameter, timevalue/3600.0))
    else:
        title = ax.set_title('Cross Section at t={0:0.2f} h'.format(timevalue/3600.0))

    title.set_position([0.5, 1.02])
    ax.set_aspect('equal')

    plt.savefig('{0}/CS_{1:05d}.png'.format(outputpath, j), dpi=300)
    print('    Outputted CS_{0:05d}.png'.format(j))

    plt.close()


def main():
    print('\nGenerating cross sections...\n')

    for j in range(outputnumber1, outputnumber2 + 1):
        index, normfactor, units = normalize(parameter)
        lists, timevalue, masslist = read_data(index, normfactor, j)

        if centering:
            lists = center_particles(lists, masslist)

        plot(lists, normfactor, units, timevalue, j)

    print('')


main()
