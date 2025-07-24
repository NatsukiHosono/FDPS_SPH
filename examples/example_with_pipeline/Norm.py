import numpy as np
import sys

#USER INPUTS ---------------------------------------------

path=str(sys.argv[1])
ncores=int(sys.argv[2])
outputnumber1=int(sys.argv[3])
outputnumber2=int(sys.argv[4])
inputparameter=str(sys.argv[5])
lower=float(sys.argv[6])
upper=float(sys.argv[7])

#---------------------------------------------------------

if inputparameter=="Density":
    parameter=9
elif inputparameter=="Energy":
    parameter=10
elif inputparameter=="Pressure":
    parameter=11
elif inputparameter=="Entropy":
    parameter=13
elif inputparameter=="Temperature":
    parameter=14

valuelist=[]
print()
for i in range(outputnumber1,outputnumber2):
    print("Reading from results."+'%05d'%i+'...')
    for j in range(ncores):
        filename=path+"results."+'%05d'%i+"_"+'%05d'%ncores+"_"'%05d'%j+".dat"
        with open(filename,"r") as file:
            file.readline()
            file.readline()
            for line in file:
                elements=line.split()
                value=float(elements[parameter])
                valuelist.append(value)

valuelist=np.array(valuelist)
lowerbound=np.percentile(valuelist,lower)
upperbound=np.percentile(valuelist,upper)

print()
print("The {} and {} percentile of {} between output numbers {} and {} are:".format(lower,upper,inputparameter,outputnumber1,outputnumber2))
print()
print("{}, {}".format(lower,upper))
print()


