import numpy as np
import scipy
import sys

try: 
    # Read the number of points from the command line if we can
    N = sys.argv[1]  
except IndexError: 
    # Otherwise use a default value
    N = 50


Phi = np.linspace(0.0, 90.0, num=N) #Lista de angulos en Phi
Theta = np.linspace(0.0, 90.0, num=N) #Lista de angulos en Theta
Radi = np.logspace(0.0, np.log10(9.0), num=N) #Lista de pasos en z de cada modelo de Cloudy
Emissivity = np.linspace(0.1, 1.0, num=N) #Lista de las emisividades de las lineas de Cloudy

NK = len(Phi)
NJ = len(Theta)
NI = len(Radi)

sumEmiss = 0.0
#La integral sobre r, th, ph
for k in range(NK):
    kneg = max(0, k - 1)
    kpos = min(NK-1, k + 1)
    DPhi = 0.5*np.radians(Phi[kpos] - Phi[kneg])
    for j in range(NJ):
        jneg = max(0, j - 1)
        jpos = min(NJ-1, j + 1)
        DMu = -0.5*(np.cos(np.radians(Theta[jpos])) - np.cos(np.radians(Theta[jneg])))
        for i in range(NI):
            ineg = max(0, i - 1)
            ipos = min(NI-1, i + 1)
            DR = 0.5*(Radi[ipos] - Radi[ineg])
            DVol =  DPhi * DMu * (Radi[i]**2) * DR
            sumEmiss += DVol

Volume = (4.*np.pi/3.) * (Radi[-1]**3 - Radi[0]**3) * 0.5 * 0.25
print "sumEmiss =", sumEmiss
print "Volume = ",  Volume
print "Relative error = ", (sumEmiss - Volume)/Volume





