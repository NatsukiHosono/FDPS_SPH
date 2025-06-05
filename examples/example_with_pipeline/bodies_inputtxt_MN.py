from __future__ import print_function, division
import os
import sys
import shutil
import stat
from math import sqrt, asin, pi
from optparse import OptionParser

def parse_args():
    usage = "usage: %prog [options] <12 parameters>"
    parser = OptionParser(usage=usage)
    parser.add_option("-c", "--config", dest="config",
                      help="Path to configuration file")

    (options, args) = parser.parse_args()

    if options.config:
        try:
            with open(options.config, 'r') as f:
                lines = [line.strip() for line in f if line.strip()]
        except IOError:
            print("Error: Cannot open configuration file.")
            sys.exit(1)
        if len(lines) != 12:
            print("Error: Configuration file must contain exactly 12 parameters.")
            sys.exit(1)
        return lines
    elif len(args) == 12:
        return args
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
with open('tar.txt', 'w') as f:
    f.write("mode = 2\n")
    f.write("UnitMass = {}\n".format(TarMass))
    f.write("UnitRadi = {}\n".format(TarRad))
    f.write("coreFracMass = {}\n".format(TarCMF))
    f.write("coreFracRadi = {}\n".format(TarCRF))
    f.write("total_number_of_particles = {}\n".format(int(Npart)))
    f.write("end_time = {}\n".format(end_time))
    f.write("damping = 1\n")
    f.write("output_interval = 100\n")
    f.write("output_directory = {}\n".format(TarOut))
    f.write("silicate_entropy = 3165.0\n")
    f.write("iron_entropy = 1500.0\n")
    f.write("impact_angle = 43\n")  # dummy parameter
    f.write("impVel = 9300\n")      # dummy parameter
    f.write("imptarMassRatio = {}\n".format(ImpMass / TarMass))  # dummy parameter

# Write impactor configuration
with open('imp.txt', 'w') as f2:
    f2.write("mode = 2\n")
    f2.write("UnitMass = {}\n".format(ImpMass))
    f2.write("UnitRadi = {}\n".format(ImpRad))
    f2.write("coreFracMass = {}\n".format(ImpCMF))
    f2.write("coreFracRadi = {}\n".format(ImpCRF))
    f2.write("total_number_of_particles = {}\n".format(int(Npart * ImpMass / TarMass)))
    f2.write("end_time = {}\n".format(end_time))
    f2.write("damping = 1\n")
    f2.write("output_interval = 100\n")
    f2.write("output_directory = {}\n".format(ImpOut))
    f2.write("silicate_entropy = 3165.0\n")
    f2.write("iron_entropy = 1500.0\n")
    f2.write("impact_angle = 43\n")  # dummy parameter
    f2.write("impVel = 9300\n")      # dummy parameter
    f2.write("imptarMassRatio = {}\n".format(ImpMass / TarMass))  # dummy parameter

# Create impactor launch script
with open('imp_relaxation.sbatch', 'w') as fb:
    fb.write("#!/bin/bash\n")
    fb.write("#SBATCH -p luna\n")
    fb.write("#SBATCH -n 100\n")
    fb.write("#SBATCH -J ImpRelax\n")
    fb.write("#SBATCH -o my_sim_%j\n")
    fb.write("#SBATCH --mem-per-cpu=5GB\n")
    fb.write("#SBATCH -t 32:00:00\n")
    fb.write("module load openmpi/4.0.3/b3\n")
    fb.write("mpirun -n 100 ./sph.out -i {}imp.txt\n".format(ImpOut))

#Create target launch script
with open('tar_relaxation.sbatch', 'w') as fa:
    fa.write("#!/bin/bash\n")
    fa.write("#SBATCH -p luna\n")
    fa.write("#SBATCH -n 100\n")
    fa.write("#SBATCH -J TarRelax\n")
    fa.write("#SBATCH -o my_sim_%j\n")
    fa.write("#SBATCH --mem-per-cpu=5GB\n")
    fa.write("#SBATCH -t 32:00:00\n")
    fa.write("module load openmpi/4.0.3/b3\n")
    fa.write("mpirun -n 100 ./sph.out -i {}tar.txt\n".format(TarOut))

# Set execute permissions and move files
os.chmod("launch_relaxation.sh", stat.S_IRWXU)
os.rename('imp.txt', os.path.join(ImpOut, 'imp.txt'))
os.rename('tar.txt', os.path.join(TarOut, 'tar.txt'))
os.rename('launch_relaxation.sh', os.path.join('..', '..', 'launch_relaxation.sh'))
