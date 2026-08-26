import numpy as np
import os, sys, shutil, stat
from math import sqrt, asin, pi

if len(sys.argv) != 13:
    print "Usage inputs: <Target N of particles> <tar M [kg]> <tar R [m]> <tar CMF> <tar CRF> <tar out path> <imp M> <imp R> <imp CMF> <imp CRF> <imp out path> <end time [s]>"
    exit(1)

Npart=np.double(sys.argv[1])
TarMass = np.double(sys.argv[2])
TarRad =np.double(sys.argv[3])
TarCMF = np.double(sys.argv[4])
TarCRF = np.double(sys.argv[5])
TarOut = sys.argv[6]
ImpMass = np.double(sys.argv[7])
ImpRad = np.double(sys.argv[8])
ImpCMF = np.double(sys.argv[9])
ImpCRF = np.double(sys.argv[10])
ImpOut = sys.argv[11]
end_time = sys.argv[12]

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

with open('imp_relaxation.sbatch','wrt') as fb:
    fb.write("#!/bin/bash\n")
    fb.write("#SBATCH -p luna\n")
    fb.write("#SBATCH -n 100\n")
    fb.write("#SBATCH -J ImpRelax\n")
    fb.write("#SBATCH -o my_sim_%j\n")
    fb.write("#SBATCH --mem-per-cpu=5GB\n")
    fb.write("#SBATCH -t 32:00:00\n")
    fb.write("module load bluehive/2.5\n")
    fb.write("module load openmpi/4.0.3/b3\n")
    fb.write("mpirun -n 100 ./sph.out -i "+ImpOut+"imp.txt\n")
fb.closed

with open('tar_relaxation.sbatch','wrt') as fa:
    fa.write("#!/bin/bash\n")
    fa.write("#SBATCH -p luna\n")
    fa.write("#SBATCH -n 100\n")
    fa.write("#SBATCH -J TarRelax\n")
    fa.write("#SBATCH -o my_sim_%j\n")
    fa.write("#SBATCH --mem-per-cpu=5GB\n")
    fa.write("#SBATCH -t 32:00:00\n")
    fa.write("module load bluehive/2.5\n")
    fa.write("module load openmpi/4.0.3/b3\n")
    fa.write("mpirun -n 100 ./sph.out -i "+TarOut+"tar.txt\n")
fa.closed

os.chmod("imp_relaxation.sbatch",stat.S_IRWXU)
os.chmod("tar_relaxation.sbatch",stat.S_IRWXU)
os.system('mv imp.txt '+ImpOut)
os.system('mv tar.txt '+TarOut)
os.system('mv imp_relaxation.sbatch ../../')
os.system('mv tar_relaxation.sbatch ../../')     
