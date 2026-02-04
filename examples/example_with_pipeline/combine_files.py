import numpy as np
import pandas as pd
import os, sys, shutil, stat


if len(sys.argv) != 5:
    print "Usage inputs: <body (tar or imp)> <body path> <output number> <ncores>"
    exit(1)

body=sys.argv[1]
path=sys.argv[2]
step=np.int(np.float(sys.argv[3]))
ncores=np.int(np.float(sys.argv[4]))

if body=="tar":

    output_path="../../input/tar.dat"

elif body=="imp":

    output_path="../../input/imp.dat"

elif body=="resume":

    output_path="../../input/resume.dat"

def concatenate(path,step,ncores,output):
    frames=[pd.read_csv(path+"results."+'%05d'%step+"_"+'%05d'%ncores+"_"'%05d'%f+".dat",skiprows=2,header=None,sep='\t') for f in range(ncores)]
    body=pd.concat(frames,ignore_index=True).sort_values(by=[0])
    with open(path+"/results."+'%05d'%step+"_"+'%05d'%ncores+"_00001.dat") as f1:
        time=f1.readline()
    f1.closed
    with open(output,'wt') as f:
        f.write(time)
        f.write(str(len(body.index))+"\n")
    f.closed
    body.to_csv(output,sep='\t',mode='a',index=False, header=False)

concatenate(path,step,ncores,output_path)
