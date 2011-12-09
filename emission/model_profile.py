import numpy as np
import scipy
import sys, os, glob
import argparse

sys.path.append("../src")       # make sure we can find claudia.py
import claudia

# Avoid verbose error messages from numpy during the reading of the Cloudy files
import warnings, argparse
warnings.filterwarnings("ignore")

# Parse command line arguments
parser = argparse.ArgumentParser(
    description="Calculate line profile of Cloudy proplyd model")

parser.add_argument(
    "--modeldir", type=str, default=".",  
    help='Path to proplyd directory model')
# parser.add_argument(
#     "--ntheta", type=int, default=100,
#     help='Number of angles theta')
parser.add_argument(
    "--nphi", type=int, default=200,
    help='Number of angles phi')
parser.add_argument(
    "--nvel", type=int, default=50,
    help='Number of velocity bins')
parser.add_argument(
    "--inc", type=float, default=0.0,
    help='Inclination in degrees')

parser.add_argument(
    "--r0", type=float, default=8.0e14,
    help='Radius of proplyd ionization front in cm')
parser.add_argument(
    "--Rmax", type=float, default=9.0,
    help='Radius of outer face of cloudy model in units of r0')

parser.add_argument(
    "--emline", type=str, default="H__1__6563A",
    help="Which emission line to use from the Cloudy model")

cmd_args = parser.parse_args()

# Set some local variables from the command line arguments
NU = cmd_args.nvel
inc_degrees = cmd_args.inc

#set de variables of the model that are defined in *.in

r0 = cmd_args.r0
Rmax = cmd_args.Rmax

thetadirs =  glob.glob(os.path.join(cmd_args.modeldir, "th??"))

# The following arrays are defined for each angle:
# Radius, Velocity, Emissivity
Radius_vectors = []
Velocity_vectors = []
Emissivity_vectors = []
for thetadir in thetadirs:
    # Read in the Cloudy model for each angle
    claudia.CloudyModel.indir = os.path.join(cmd_args.modeldir, thetadir)
    claudia.CloudyModel.outdir = claudia.CloudyModel.indir
    # find the last iteration
    modelname = glob.glob(os.path.join(claudia.CloudyModel.indir, "*.in"))[-1]
    m = claudia.CloudyModel(modelname) 
    Radius = m.ovr.depth #Lista de pasos en z de cada modelo de Cloudy
    R = Rmax - Radius/r0            # array of dimensionless radius
    i2 = len(R[R>=1.0])             # find index where R = 1

    # Calculate physical variables as function of radius for this theta
    Density = 10**(m.ovr.hden)
    Emissivity = 10**(m.em[cmd_args.emline])
    sound_speed = np.sqrt(3./5.)*m.pre.cadwind_kms # isothermal sound speed
    n0 = Density[i2]
    Velocity = sound_speed[i2] * n0 / (R**2 * Density) 

    # Add them to the lists of vectors 
    Radius_vectors.append(Radius)
    Velocity_vectors.append(Velocity)
    Emissivity_vectors.append(Emissivity)
    

#The model

    
cosi = np.cos(np.radians(inc_degrees))
sini = np.sin(np.radians(inc_degrees))

# All arrays begin with capital letters
phimin, phimax = 0.0, np.radians(360.0)
thetamin, thetamax = 0.0, np.radians(90.0)

Phi = np.linspace(phimin, phimax, num=cmd_args.nphi) #Lista de angulos en Phi
# For the moment, we use the existing thetas with no interpolation
Theta = np.array([float(s[2:]) for s in thetadirs])
# Theta = np.linspace(thetamin, thetamax, num=cmd_args.ntheta) #Lista de angulos en Theta
NK = len(Phi)
NJ = len(Theta)


# Define the velocity bins
umax = max(abs(Velocity_vectors[0]))
umin = -umax
du = (umax - umin)/NU
Perfil = np.zeros(NU)

sumEmiss = 0.0

#Line profile

for k in range(NK):
    cphi = np.cos(Phi[k])
    sphi = np.sin(Phi[k])
    kneg = max(0, k - 1)
    kpos = min(NK-1, k + 1)
    dphi = 0.5*(Phi[kpos] - Phi[kneg])
    for j in range(NJ):
        # for each angle
        Radius = Radius_vectors[j]
        Velocity = Velocity_vectors[j]
        Emissivity = Emissivity_vectors[j]

        ctheta = np.cos(Theta[j])
        stheta = np.sin(Theta[j])
        jneg = max(0, j - 1)
        jpos = min(NJ-1, j + 1)
        dmu = -0.5*(np.cos(Theta[jpos]) - np.cos(Theta[jneg]))
        NI = len(Radius)
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
PerfilU = np.linspace(umin, umax, NU)
savefile = "%(modelname)s-perfil-ntheta%(ntheta)i-nphi%(nphi)i-nvel%(nvel)i-inc%(inc)i-%(emline)s.dat" % (vars(cmd_args))
np.savetxt(savefile, (PerfilU, Perfil))

Volume = (4.*np.pi/3.) * (Radius[-1]**3 - Radius[0]**3) * 0.5
print "sumEmiss =", sumEmiss
print "Volume = ",  Volume
print "Relative error = ", (sumEmiss - Volume)/Volume

print "Sum of line profile (should be same as sumEmiss): ", Perfil.sum()
print Perfil


