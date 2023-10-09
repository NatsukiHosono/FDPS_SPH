import numpy as np
import pandas as pd
import os, sys, shutil, stat


if len(sys.argv) != 6:
    print "Usage inputs: <target path> <impactor path> <ncores> <tar time> <impactor time>"
    exit(1)

paths=[sys.argv[1],sys.argv[2]]
ncores=np.int(np.float(sys.argv[3]))
steps=[np.int(np.float(sys.argv[4]))/100,np.int(np.float(sys.argv[5]))/100]
output_paths=["../../input/tar.dat","../../input/imp.dat"]

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

for i in range(len(paths)):
     concatenate(paths[i],steps[i],ncores,output_paths[i])

