import os, sys, shutil, stat
from math import sqrt, asin, pi
import matplotlib.pyplot as plt
import numpy as np
from minimaldefinitions import ReverseTime

if len(sys.argv) != 6:
    print "Usage inputs: <b> <vimp [v_esc]> <initial separation [R_tar]> <end_time [s]> <output_path>" 
    exit(1)

b = np.double(sys.argv[1])
v_esc_multiple=np.double(sys.argv[2])
separation = np.double(sys.argv[3])
end_time = np.double(sys.argv[4])
output_path = sys.argv[5]

name = "5_new"
target_path = "../../input/tar.dat"
impactor_path = "../../input/imp.dat"
dt = 5
rev_output_path = output_path+"{}_{}_reversed_outputs".format(name, b)

if os.path.exists(output_path):
    shutil.rmtree(output_path)
os.mkdir(output_path)

rt = ReverseTime(
    target_file_path=target_path,
    impactor_file_path=impactor_path,
    impact_parameter=b,
    dt=-dt,
    center_target=True,
    v_esc_multiple=v_esc_multiple
)

v_esc = sqrt(2.0 * rt.G * (rt.target_mass + rt.impactor_mass) / (rt.radius_target + rt.radius_impactor))
v_imp = v_esc_multiple * v_esc

loop = 0
plotted = 0
while rt.distance < separation * rt.radius_target:
    rt.reverse()
    if loop % 5 == 0:
#        fig = rt.plot_current_position()
#        plt.savefig(output_path + "/{}.png".format(plotted), format='png')
        plotted += 1
    loop += 1

# fig = rt.plot_current_position()
# plt.show()

print(
    "INITIAL SETUP:\n"
    "TARGET PATH: {}\n"
    "IMPACTOR PATH: {}\n"
    "b: {} ({} degrees)\n"
    "dt: -{} s\n"
    "IMPACTOR-TARGET DISTANCE: {} Rad Tar ({} m)\n"
    "IMPACTOR-TARGET X DISTANCE: {} Rad Tar ({} m)\n"
    "IMPACTOR-TARGET Y DISTANCE: {} Rad Tar ({} m)\n"
    "IMPACTOR-TARGET Z DISTANCE: {} Rad Tar ({} m)\n"
	    "INITIAL TARGET COORDS: {}\n"
    "INITIAL IMPACTOR COORDS: {}\n"
    "ESCAPE VELOCITY: {}\n"
    "IMPACT VELOCITY: {}\n" 
    "TARGET RADIUS : {}\n"
    "IMPACTOR RADIUS: {}\n"
    "TARGET MASS: {}\n"
    "IMPACTOR MASS: {}\n"
    "INITIAL_TARGET_X_VELOCITY: {} v_esc\n"
    "INITIAL_TARGET_Y_VELOCITY: {} v_esc\n"
    "INITIAL_TARGET_Z_VELOCITY: {} v_esc\n"
    "INITIAL_IMPACTOR_X_VELOCITY: {} v_esc\n"
    "INITIAL_IMPACTOR_Y_VELOCITY: {} v_esc\n"
    "INITIAL_IMPACTOR_Z_VELOCITY: {} v_esc\n".format(
        target_path, impactor_path,
        b, asin(b) * (180 / pi),
        dt,
        rt.distance / rt.radius_target,
        rt.x_distance / rt.radius_target, rt.y_distance / rt.radius_target, rt.z_distance / rt.radius_target,
        rt.distance, rt.x_distance, rt.y_distance, rt.z_distance,
        rt.com_target, rt.com_impactor, v_esc, v_imp, rt.radius_target, rt.radius_impactor, rt.target_mass, rt.impactor_mass,
        rt.v_target_x / v_esc, rt.v_target_y / v_esc, rt.v_target_z / v_esc,
        rt.v_impactor_x / v_esc, rt.v_impactor_y / v_esc, rt.v_impactor_z / v_esc,
    )
)
'''
animate(
    start_time=0,
    end_time=plotted - 1,
    interval=1,
    path=output_path,
    fps=10,
    filename="{}_reverse_time.mp4".format(name),
    reverse=True
)
'''
rt.plot_velocity_history()

with open('gi.txt','wt') as f:
    f.write("mode = 1\n")
    f.write("end_time = "+str(end_time)+"\n")
    f.write("damping = 1\n")
    f.write("output_interval = 100\n")
    f.write("output_directory = "+output_path+"\n")
    f.write("impact_angle = "+str(asin(b)*180/pi)+"\n")
    f.write("L_init_vs_L_em = 0.1\n")
    f.write("silicate_entropy = 3000.0\n")
    f.write("iron_entropy = 1750.0\n")
    f.write("silicate_grid_size = 120\n")
    f.write("iron_grid_size = 120\n")
    f.write("v_imp_x_v_esc = "+str(rt.v_impactor_x/v_esc)+"\n")
    f.write("v_imp_y_v_esc = "+str(rt.v_impactor_y/v_esc)+"\n")
    f.write("v_imp_z_v_esc = "+str(rt.v_impactor_z/v_esc)+"\n")
    f.write("v_tar_x_v_esc = "+str(rt.v_target_x/v_esc)+"\n")
    f.write("v_tar_y_v_esc = "+str(rt.v_target_y/v_esc)+"\n")
    f.write("v_tar_z_v_esc = "+str(rt.v_target_z/v_esc)+"\n")
    f.write("imp_x_init = "+str(rt.com_impactor[0])+"\n")
    f.write("imp_y_init = "+str(rt.com_impactor[1])+"\n")
    f.write("imp_z_init = "+str(rt.com_impactor[2])+"\n")
    f.write("tar_x_init = 0.0\n")
    f.write("tar_y_init = 0.0\n")
    f.write("tar_z_init = 0.0\n")
f.closed

with open('launch_gi.sh','wrt') as fb:
    fb.write("#!/bin/bash\n")
    fb.write("#SBATCH -p luna\n")
    fb.write("#SBATCH -n 100\n")
    fb.write("#SBATCH -J TarImpRelax\n")
    fb.write("#SBATCH -t 32:00:00\n")
    fb.write("module load openmpi/4.0.3/b3\n")
    fb.write("mpirun -np 100 ./sph.out -i "+output_path+"gi.txt\n")
fb.closed

os.chmod("launch_gi.sh",stat.S_IRWXU)
os.system('mv gi.txt '+output_path)
os.system('mv launch_gi.sh ../../')

