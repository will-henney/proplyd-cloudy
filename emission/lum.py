import numpy as np
import scipy
import sys

try: 
    # Read the number of points from the command line if we can
    N = sys.argv[1]  
except IndexError: 
    # Otherwise use a default value
    N = 50
try: 
    # Read the number of bins from the command line if we can
    NU = sys.argv[2]  
except IndexError: 
    # Otherwise use a default value
    NU = 50
try: 
    # Read the projection angle if we can
    inclination = sys.argv[3]  
except IndexError:
    # Otherwise use a default value
    inclination = 0


# All arrays begin with capital letters
Phi = np.linspace(0.0, 90.0, num=N) #Lista de angulos en Phi
Theta = np.linspace(0.0, 90.0, num=N) #Lista de angulos en Theta
Radius = np.logspace(0.0, np.log10(9.0), num=N) #Lista de pasos en z de cada modelo de Cloudy
Emissivity = np.linspace(0.1, 1.0, num=N) #Lista de las emisividades de las lineas de Cloudy
Velocity = 20*(np.ones(N))

# Define the velocity bins
umax = max(Velocity)
umin = -min(Velocity)
DU = (umax - umin)/NU
Perfil=np.zeros(NU)

NK = len(Ph)
NJ = len(Th)
NI = len(Radi)

sumEmiss = 0.0
#La integral sobre r, th, ph
for k in range(NK):
    kneg = max(0, k - 1)
    kpos = min(NK-1, k + 1)
    DPhi = 0.5*np.radians(Ph[kpos] - Ph[kneg])
    for j in range(NJ):
        jneg = max(0, j - 1)
        jpos = min(NJ-1, j + 1)
        DMu = -0.5*(np.cos(np.radians(Th[jpos])) - np.cos(np.radians(Th[jneg])))
        for i in range(NI):
            ineg = max(0, i - 1)
            ipos = min(NI-1, i + 1)
            DR = 0.5*(Radi[ipos] - Radi[ineg])
            DVol =  DPhi * DMu * (Radi[i]**2) * DR
            sumEmiss += DVol
            u = -Velocity[i]*((np.sin(np.radians(Th))+np.cos(np.radians(i)))*(np.cos(np.radians(Th))+np.sin(np.radians(Th))*np.cos(np.radians(Ph)))+(np.sin(np.radians(Th))*np.sin(np.radians(Ph))))
            x = (u-umin)/DU
            I = int(x)

Volume = (4.*np.pi/3.) * (Radi[-1]**3 - Radi[0]**3) * 0.5 * 0.25
print "sumEmiss =", sumEmiss
print "Volume = ",  Volume
print "Relative error = ", (sumEmiss - Volume)/Volume





