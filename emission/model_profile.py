import numpy as np
import scipy
import sys, os, glob
import argparse

sys.path.append("../../../src")       # make sure we can find claudia.py
import claudia

# Avoid verbose error messages from numpy during the reading of the Cloudy files
import warnings, argparse
warnings.filterwarnings("ignore")

# Parse command line arguments
parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="""
Calculate 1-d line profile of Cloudy proplyd model. 
If the --window option is used, 
also calculate the 3-d position-position-velocity cube.
""")

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
    "--window", type=float, nargs=4, default=argparse.SUPPRESS,
    metavar=("X1", "Y1", "X2", "Y2"),
    help='If present, bounding box of view window (in units of r0) for the position-position-velocity cube')
parser.add_argument(
    "--pixel-size", type=float, default=0.1, 
    help='Size of pixels (in units of r0) for the PPV cube')

cmd_args = parser.parse_args()

# Set some local variables from the command line arguments
NU = cmd_args.nvel
inc_degrees = cmd_args.inc

#set de variables of the model that are defined in *.in

r0 = cmd_args.r0
Rmax = cmd_args.Rmax

thetadirs =  [
    os.path.basename(p) for p in 
    glob.glob(os.path.join(cmd_args.modeldir, "th??"))
    ]

# The following arrays are defined for each angle:
# Radius, Velocity, Emissivity
Radius_vectors = []
Velocity_vectors = []
Emissivities_vectors = []
# List of names of emission lines
emlines = None
for thetadir in thetadirs:
    # Read in the Cloudy model for each angle
    claudia.CloudyModel.indir = os.path.join(cmd_args.modeldir, thetadir)
    claudia.CloudyModel.outdir = claudia.CloudyModel.indir
    # find the last iteration
    modelinputpath = glob.glob(os.path.join(claudia.CloudyModel.indir, "*.in"))[-1]
    # extract only the basename of the file, without .in extension
    modelname = os.path.splitext(os.path.basename(modelinputpath))[0]
    print "Indir: ", claudia.CloudyModel.indir
    print "Modelname: ", modelname
    m = claudia.CloudyModel(modelname) 
    Radius = Rmax*r0 - m.ovr.depth # physical radia from proplyd center
    R = Radius/r0            # array of dimensionless radius
    i2 = len(R[R>=1.0])             # find index where R = 1

    # Find the list of emission lines if necesssary
    if emlines is None:
        emlines = m.em.dtype.names[1:] # miss off "depth"

    # Calculate physical variables as function of radius for this theta
    Density = 10**(m.ovr.hden)
    Emissivities = dict()
    for emline in emlines:
        Emissivities[emline] = 10**(m.em[emline])
    sound_speed = np.sqrt(3./5.)*m.pre.cadwind_kms # isothermal sound speed
    n0 = Density[i2]
    Velocity = sound_speed[i2] * n0 / (R**2 * Density) 

    # Add them to the lists of vectors 
    Radius_vectors.append(Radius)
    Velocity_vectors.append(Velocity)
    Emissivities_vectors.append(Emissivities)
    

#The model

    
cosi = np.cos(np.radians(inc_degrees))
sini = np.sin(np.radians(inc_degrees))

# All arrays begin with capital letters
phimin, phimax = 0.0, np.radians(360.0)
thetamin, thetamax = 0.0, np.radians(90.0)

Phi = np.linspace(phimin, phimax, num=cmd_args.nphi) #Lista de angulos en Phi
# For the moment, we use the existing thetas with no interpolation
Theta = np.array([np.radians(float(s[2:])) for s in thetadirs])
print "Theta: ", Theta
# Theta = np.linspace(thetamin, thetamax, num=cmd_args.ntheta) #Lista de angulos en Theta
NK = len(Phi)
NJ = len(Theta)


# Define the velocity bins
umax = max(abs(Velocity_vectors[0]))
umin = -umax
du = (umax - umin)/NU
Perfiles = dict()
Fluxes = dict()
for emline in emlines:
    # Line profile as function of velocity
    Perfiles[emline] = np.zeros(NU)
    # Total line flux (integral of line profile)
    Fluxes[emline] = 0.0

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
        Emissivities = Emissivities_vectors[j]

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
            u = -Velocity[i]*(sini*stheta*cphi + cosi*ctheta)
            iu = int((u-umin)/du)
            # assert iu >= 0 and iu < NU, "Index (%i) out of bounds [%i:%i] of velocity array. u = %.3g, umin, umax = %.3g, %.3g" % (iu, 0, NU-1, u, umin, umax)
            # for each emission lines
            for emline in emlines:
                # add in to line profile
                if iu >= 0 and iu < NU: 
                    Perfiles[emline][iu] += dvol * Emissivities[emline][i]
                # and add in to total flux
                Fluxes[emline] += dvol * Emissivities[emline][i]


PerfilU = np.linspace(umin, umax, NU)
savefile = "model_profile-nphi%(nphi)i-nvel%(nvel)i-inc%(inc)i.dat" % (vars(cmd_args))

import matplotlib.mlab as mlab
datatable = np.rec.fromarrays([PerfilU] + Perfiles.values(), 
                              names=','.join(['U'] + Perfiles.keys()))
mlab.rec2csv(datatable, savefile, delimiter="\t")


hbeta = Fluxes["H__1__4861A"]
print "H beta line flux: "
print "Line fluxes relative to H beta = 100:"

def sortkey(cloudy_line):
    """
    Key for sorting the lines by wavelength
    """
    cloudy_line = cloudy_line.replace('_', ' ')
    # strip off first 4 characters and then only keep the first word
    wav_s = cloudy_line[4:].split()[0]
    unit = wav_s[-1]
    wav = float(wav_s[:-1])
    if unit == 'A':
        wav *= 1.e-8
    elif unit == 'm':
        wav *= 1.e-4
    else:
        raise ValueError, "%s is neither A nor m" % (unit)
    return wav

for emline in sorted(emlines, key=sortkey):
    flux = Fluxes[emline]
    em = emline.replace('_', ' ')
    print '|'.join(['', em[:-5], em[-5:], "%g" % (100.0*flux/hbeta), ''])



