import numpy as np

Phi = [0,10,20,30,40,50,60,70,80,90] #Lista de angulos en Phi
Theta = [0,10,20,30,40,50,60,70,80,90] #Lista de angulos en Theta
Radi = [0,1,2,3,4,5,6,7,8,9] #Lista de pasos en z de cada modelo de Cloudy
Emissivity = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,10.] #Lista de las emisividades de las lineas de Cloudy

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

# #La integral sobre Theta

# for j in range(NJ):
#     sumTheta = sumTheta + DTheta
#     print "sumTheta =", sumTheta

# #La integral sobre Phi

# for k in range(NK):
#     kneg = max(0, k - 1)
#     kpos = min(NK-1, k + 1)
#     DPhi = 0.5*(Phi[kpos] - Phi[kneg])
#     sumPhi = sumPhi + DPhi
#     print "sumPhi =", sumPhi


