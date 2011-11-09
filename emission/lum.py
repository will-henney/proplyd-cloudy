import numpy as np
import scipy

Phi = scipy.linspace(0, 90, num=100) #Lista de angulos en Phi
Theta = scipy.linspace(0, 90, num=100) #Lista de angulos en Theta
Radi = scipy.linspace(0, 9, num=1000) #Lista de pasos en z de cada modelo de Cloudy
Emissivity = scipy.linspace(0.1, 1, num=100) #Lista de las emisividades de las lineas de Cloudy

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

print "sumEmiss =", sumEmiss
print "Volume = ",  (4.*np.pi/3.) * (Radi[-1]**3 - Radi[0]**3) * 0.5 * 0.25




