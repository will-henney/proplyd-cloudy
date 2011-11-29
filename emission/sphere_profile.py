import numpy as np
import scipy
import sys, os
import argparse

sys.path.append("../src")       # make sure we can find claudia.py
import claudia


# Parse command line arguments
parser = argparse.ArgumentParser(
    description="Calculate line profile of spherical Cloudy model")

parser.add_argument(
    "modelpath", type=str,
    help='Path to Cloudy model (assumes input and output in same directory)')
parser.add_argument(
    "--ntheta", type=int, default=50,
    help='Number of angles theta')
parser.add_argument(
    "--nphi", type=int, default=200,
    help='Number of angles phi')
parser.add_argument(
    "--nvel", type=int, default=50,
    help='Number of velocity bins')
parser.add_argument(
    "--inc", type=float, default=0.0,
    help='Inclination in degrees')

cmd_args = parser.parse_args()

# Set some local variables from the command line arguments
NU = cmd_args.nvel
inc_degrees = cmd_args.inc
modeldir, modelname = os.path.split(cmd_args.modelpath)

# Read in Cloudy model
if modeldir != '':
    claudia.CloudyModel.indir = modeldir
    claudia.CloudyModel.outdir = modeldir
else:
    claudia.CloudyModel.indir = '.'
    claudia.CloudyModel.outdir = '.'
m = claudia.CloudyModel(modelname) 
    
cosi = np.cos(np.radians(inc_degrees))
sini = np.sin(np.radians(inc_degrees))

# All arrays begin with capital letters
phimin, phimax = 0.0, np.radians(360.0)
thetamin, thetamax = 0.0, np.radians(90.0)

Phi = np.linspace(phimin, phimax, num=cmd_args.nphi) #Lista de angulos en Phi
Theta = np.linspace(thetamin, thetamax, num=cmd_args.ntheta) #Lista de angulos en Theta
Radius = m.ovr.depth #Lista de pasos en z de cada modelo de Cloudy
NK = len(Phi)
NJ = len(Theta)
NI = len(Radius)

Velocity = m.pre.cadwind_kms
Emissivity = m.str

# Define the velocity bins
umax = max(abs(Velocity))
umin = -umax
du = (umax - umin)/NU
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
            sumEmiss += dvol * Emissivity[i]
            u = -Velocity[i]*(sini*stheta*cphi + cosi*ctheta)
            x = (u-umin)/du
            iu = int(x)
            assert iu >= 0 and iu < NU, "Index (%i) out of bounds [%i:%i] of velocity array" % (iu, 0, NU-1)
            Perfil[iu] += dvol * Emissivity[i]

Volume = (4.*np.pi/3.) * (Radius[-1]**3 - Radius[0]**3) * 0.5
print "sumEmiss =", sumEmiss
print "Volume = ",  Volume
print "Relative error = ", (sumEmiss - Volume)/Volume

print "Sum of line profile (should be same as sumEmiss): ", Perfil.sum()
print Perfil

PerfilU = np.linspace(umin, umax, NU)
savefile = "%(modelname)s-perfil-NU%(NU)i-inc%(inc_degrees)i.dat" % (locals())
# savefile = "%(modelname)s-perfil-N%(N)i-NU%(NU)i-inc%(inc_degrees)i.dat" % (locals())
np.savetxt(savefile, (PerfilU, Perfil))
