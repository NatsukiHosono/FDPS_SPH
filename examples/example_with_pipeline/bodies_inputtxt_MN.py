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

for p in [TarOut, ImpOut]:   
    if os.path.exists(p):
        shutil.rmtree(p)
    os.mkdir(p)

with open('tar.txt','wt') as f:
    f.write("mode = 2\n")
    f.write("UnitMass = "+str(TarMass)+"\n")
    f.write("UnitRadi = "+str(TarRad)+"\n")
    f.write("coreFracMass = "+str(TarCMF)+"\n")
    f.write("coreFracRadi = "+str(TarCRF)+"\n")
    f.write("total_number_of_particles ="+str(int(Npart))+"\n")
    f.write("end_time = "+str(end_time)+"\n")
    f.write("damping = 1\n")
    f.write("output_interval = 100\n")
    f.write("output_directory = "+TarOut+"\n")
    f.write("silicate_entropy = 3165.0\n")
    f.write("iron_entropy = 1500.0\n")
    f.write("impact_angle = 43\n") #dummy parameter
    f.write("impVel = 9300\n") #dummy parameter
    f.write("imptarMassRatio = "+str(ImpMass/TarMass)+"\n") #dummy parameter
f.closed

with open('imp.txt','wt') as f2:
    f2.write("mode = 2\n")
    f2.write("UnitMass = "+str(ImpMass)+"\n")
    f2.write("UnitRadi = "+str(ImpRad)+"\n")
    f2.write("coreFracMass = "+str(ImpCMF)+"\n")
    f2.write("coreFracRadi = "+str(ImpCRF)+"\n")
    f2.write("total_number_of_particles = "+str(int(Npart*ImpMass/TarMass))+"\n")
    f2.write("end_time = "+str(end_time)+"\n")
    f2.write("damping = 1\n")
    f2.write("output_interval = 100\n")
    f2.write("output_directory = "+ImpOut+"\n")
    f2.write("silicate_entropy = 3165.0\n")
    f2.write("iron_entropy = 1500.0\n")
    f2.write("impact_angle = 43\n") #dummy parameter
    f2.write("impVel = 9300 \n") #dummy parameter
    f2.write("imptarMassRatio = "+str(ImpMass/TarMass)+"\n") #dummy parameter
f2.closed

with open('launch_relaxation.sh','wrt') as fb:
    fb.write("#!/bin/bash\n")
    fb.write("#SBATCH -p luna\n")
    fb.write("#SBATCH -n 100\n")
    fb.write("#SBATCH -J TarImpRelax\n")
    fb.write("#SBATCH -t 32:00:00 -o out.%a.txt -a 1-2\n")
    fb.write("module load openmpi/4.0.3/b3\n")
    fb.write("mpirun -n 100 ./sph.out -i "+TarOut+"tar.txt\n")
    fb.write("mpirun -n 100 ./sph.out -i "+ImpOut+"imp.txt\n")
fb.closed

os.chmod("launch_relaxation.sh",stat.S_IRWXU)
os.system('mv imp.txt '+ImpOut)
os.system('mv tar.txt '+TarOut)
os.system('mv launch_relaxation.sh ../../')            
