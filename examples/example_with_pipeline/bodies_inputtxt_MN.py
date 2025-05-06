import numpy as np
import os
import sys
import shutil
import stat
import argparse
from math import sqrt, asin, pi

def parse_args():
    parser = argparse.ArgumentParser(description='Process simulation parameters.')
    parser.add_argument('--config', type=str, help='Path to configuration file')
    parser.add_argument('params', nargs='*', help='Manual input parameters')

    args = parser.parse_args()

    if args.config:
        with open(args.config, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]
        if len(lines) != 12:
            print("Error: Configuration file must contain exactly 12 parameters.")
            sys.exit(1)
        return lines
    elif len(args.params) == 12:
        return args.params
    else:
        print("Error: Provide either a configuration file with 12 parameters or 12 manual parameters.")
        sys.exit(1)

# Parse arguments
params = parse_args()

# Assign parameters
Npart = float(params[0])
TarMass = float(params[1])
TarRad = float(params[2])
TarCMF = float(params[3])
TarCRF = float(params[4])
TarOut = params[5]
ImpMass = float(params[6])
ImpRad = float(params[7])
ImpCMF = float(params[8])
ImpCRF = float(params[9])
ImpOut = params[10]
end_time = params[11]

# Create output directories
for p in [TarOut, ImpOut]:   
    if os.path.exists(p):
        shutil.rmtree(p)
    os.mkdir(p)

# Write target configuration
with open('tar.txt', 'wt') as f:
    f.write("mode = 2\n")
    f.write(f"UnitMass = {TarMass}\n")
    f.write(f"UnitRadi = {TarRad}\n")
    f.write(f"coreFracMass = {TarCMF}\n")
    f.write(f"coreFracRadi = {TarCRF}\n")
    f.write(f"total_number_of_particles = {int(Npart)}\n")
    f.write(f"end_time = {end_time}\n")
    f.write("damping = 1\n")
    f.write("output_interval = 100\n")
    f.write(f"output_directory = {TarOut}\n")
    f.write("silicate_entropy = 3165.0\n")
    f.write("iron_entropy = 1500.0\n")
    f.write("impact_angle = 43\n")  # dummy parameter
    f.write("impVel = 9300\n")      # dummy parameter
    f.write(f"imptarMassRatio = {ImpMass / TarMass}\n")  # dummy parameter

# Write impactor configuration
with open('imp.txt', 'wt') as f2:
    f2.write("mode = 2\n")
    f2.write(f"UnitMass = {ImpMass}\n")
    f2.write(f"UnitRadi = {ImpRad}\n")
    f2.write(f"coreFracMass = {ImpCMF}\n")
    f2.write(f"coreFracRadi = {ImpCRF}\n")
    f2.write(f"total_number_of_particles = {int(Npart * ImpMass / TarMass)}\n")
    f2.write(f"end_time = {end_time}\n")
    f2.write("damping = 1\n")
    f2.write("output_interval = 100\n")
    f2.write(f"output_directory = {ImpOut}\n")
    f2.write("silicate_entropy = 3165.0\n")
    f2.write("iron_entropy = 1500.0\n")
    f2.write("impact_angle = 43\n")  # dummy parameter
    f2.write("impVel = 9300\n")      # dummy parameter
    f2.write(f"imptarMassRatio = {ImpMass / TarMass}\n")  # dummy parameter

# Create launch script
with open('launch_relaxation.sh', 'wt') as fb:
    fb.write("#!/bin/bash\n")
    fb.write("#SBATCH -p luna\n")
    fb.write("#SBATCH -n 100\n")
    fb.write("#SBATCH -J TarImpRelax\n")
    fb.write("#SBATCH -t 32:00:00 -o out.%a.txt -a 1-2\n")
    fb.write("module load openmpi/4.0.3/b3\n")
    fb.write(f"mpirun -n 100 ./sph.out -i {TarOut}tar.txt\n")
    fb.write(f"mpirun -n 100 ./sph.out -i {ImpOut}imp.txt\n")

# Set execute permissions and move files
os.chmod("launch_relaxation.sh", stat.S_IRWXU)
os.system(f'mv imp.txt {ImpOut}')
os.system(f'mv tar.txt {TarOut}')
os.system('mv launch_relaxation.sh ../../')
