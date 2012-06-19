import numpy as np
import pyfits
import matplotlib.pyplot as plt
import argparse
import sys

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="""
Plot line profiles from 1D FITS files that were created by extract_aperture.py
    """)

parser.add_argument(
    "lines", type=str, nargs="+",  
    help='Cloudy IDs of emission lines to plot (with _ instead of space)')

parser.add_argument(
    "--raw", action="store_true",
    help='Use the raw fluxes (default is to normalize to peak)')

cmd_args = parser.parse_args()

for line in cmd_args.lines:
    fitsfile = "profile-{}.fits".format(line)
    hdu, = pyfits.open(fitsfile)

    k0, u0, du, nu = [hdu.header[kwd] for kwd in "CRPIX1", "CRVAL1", "CDELT1", "NAXIS1"]
    U = u0 + du*(np.arange(nu) - k0) # velocity array
    F = hdu.data                     # flux array
    if not cmd_args.raw:
        F /= F.max()
    plt.plot(U, F, label=line)

figfile = "plot_profile.png"
plt.xlabel("Velocity, km/s")
if cmd_args.raw:
    plt.ylabel("Raw flux")
else:
    plt.ylabel("Normalized flux")

plt.legend()
plt.savefig(figfile)
