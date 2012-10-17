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
        "--modeldir", type=str, default=".",  
        help='Path to proplyd directory model')
parser.add_argument(
        "--rebin", type=str, default="linear-101",
        help='Load rebinned model from "rebin-REBIN"')  
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
dmu = Mu[1] - Mu[0]             # assume uniform spacing in mu
R = read_fits("R")[0,:]
i2 = len(R[R<1.0])             # find index where R = 1
Density2D = read_fits("ovr-hden")
validmask = Density2D != -1.0
Density2D = 10**Density2D
sound_speed2D = np.sqrt(3./5.)*read_fits("pre-cadwind_kms")
n0 = Density2D[:, i2]
c0 = sound_speed2D[:, i2] * 1.0e5 # sound speed in cm/s
mp = 1.67262158e-24                # proton mass in g

loss = 2*np.pi*1.3*mp*dmu*np.sum(n0 * c0) * cmd_args.r0**2

print "Mass loss rate: ", loss, "g/s"
loss *= 3.15576e7 / 1.989e33
print "Mass loss rate: ", loss, "Msun/yr"


