import pandas as pd
import numpy as np
import scipy.constants
import os as os

filename = '2ME_b07_cold'
pathinput = '../../output/Analysis/' + filename +'/'

gi=pd.read_csv(pathinput+filename+'_identified_new.dat',sep='\t')

disk=gi.loc[gi['tag2']==1]
escaping=gi.loc[gi['tag2']==2]
fullplanet = gi.loc[gi['tag2']==0]

m = fullplanet['m'].sum()
md = disk['m'].sum()
ME =  5.9722e24

omega = ((planet['jx'].sum()**2+planet['jy'].sum()**2+planet['jz'].sum()**2))**0.5/(planet['x']**2+planet['y']**2).sum()

print("Masses of planet and disk: [ME]", m/ME, md/ME)
print("Omega, rotation period [h]: ", omega, 2*np.pi/omega/(3600))