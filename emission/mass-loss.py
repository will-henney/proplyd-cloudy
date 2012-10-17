import pyfits
import glob
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="""
Make the integration to find the mass loss rate in the proplyd head
""")

parser.add_argument(
        "--rebin", type=str, default="linear-101",
        help='Load rebinned model from "rebin-REBIN"')  
    help='Last part of cube file names')
parser.add_argument(
        "--r0", type=float, default=8.0e14,
        help='Radius of proplyd ionization front in cm')

cmd_args = parser.parse_args()

def read_fits(varid):
    """
       Read variable VARID from fits file of interpolated model structure
       """
    fitsname = os.path.join(
           cmd_args.modeldir,
           "rebin-%s" % (cmd_args.rebin),
           "%s.fits" % (varid)
           )
    return pyfits.open(fitsname)[0].data


Mu = read_fits("mu")[:,0]
dmu = 
R = read_fits("R")[0,:]
i2 = len(R[R<1.0])             # find index where R = 1
Density2D = read_fits("ovr-hden")
validmask = Density2D != -1.0
Density2D = 10**Density2D
sound_speed2D = np.sqrt(3./5.)*read_fits("pre-cadwind_kms")
n0 = Density2D[:, i2:i2+1]
c0 = sound_speed2D[:, i2:i2+1]
R0 = R[:, i2:i2+1]
mp = 1.67262158e-24                # proton mass in g

for mu in Mu:
    dloss = n0 * c0 * (R0*r0)^2
    loss += 2*np.pi*mp*dloss*dmu
