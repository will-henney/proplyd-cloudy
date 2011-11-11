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
    inc = sys.argv[3]  
except IndexError:
    # Otherwise use a default value
    inc = 0
    
cosi = np.cos(inc)
sini = np.sin(inc)

# All arrays begin with capital letters
phimin, phimax = 0.0, np.radians(360.0)
thetamin, thetamax = 0.0, np.radians(90.0)

Phi = np.linspace(phimin, phimax, num=4*N) #Lista de angulos en Phi
Theta = np.linspace(thetamin, thetamax, num=N) #Lista de angulos en Theta
Radius = np.logspace(0.0, np.log10(9.0), num=N) #Lista de pasos en z de cada modelo de Cloudy
NK = len(Phi)
NJ = len(Theta)
NI = len(Radius)

Velocity = 20*(np.ones(NI))
Emisivity = np.ones(NI)

# Define the velocity bins
umax = max(abs(Velocity))
umin = -umax
DU = (umax - umin)/NU
Perfil=np.zeros(NU)

sumEmiss = 0.0
#La integral sobre r, th, ph
for k in range(NK):
    cphi = np.cos(Phi[k])
    sphi = np.sin(Phi[k])
    kneg = max(0, k - 1)
    kpos = min(NK-1, k + 1)
    dphi = 0.5*(Phi[kpos] - Phi[kneg])
    for j in range(NJ):
        ctheta = np.cos(Theta[j])
        stheta = np.sin(Theta[j])
        jneg = max(0, j - 1)
        jpos = min(NJ-1, j + 1)
        dmu = -0.5*(np.cos(Theta[jpos]) - np.cos(Theta[jneg]))
        for i in range(NI):
            ineg = max(0, i - 1)
            ipos = min(NI-1, i + 1)
            dr = 0.5*(Radius[ipos] - Radius[ineg])
            dvol =  dphi * dmu * (Radius[i]**2) * dr
            sumEmiss += dvol
            u = -Velocity[i]*(sini*stheta*cphi + cosi*ctheta)
            x = (u-umin)/DU
            I = int(x)

Volume = (4.*np.pi/3.) * (Radius[-1]**3 - Radius[0]**3) * 0.5
print "sumEmiss =", sumEmiss
print "Volume = ",  Volume
print "Relative error = ", (sumEmiss - Volume)/Volume





