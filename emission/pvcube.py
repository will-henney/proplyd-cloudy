"""
Swap the axes of the PPV cubes to be VPP so that SAOimage can show them as PV diagrams
"""

import pyfits
import glob
import numpy as np
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="""
Swap the axes of the PPV cubes to be VPP so that SAOimage can show them as PV diagrams
    """)

parser.add_argument(
    "--suffix", type=str, default="nphi400-nvel200-inc60",  
    help='Last part of cube file names')

cmd_args = parser.parse_args()


files = glob.glob("cube-*-{}.fits".format(cmd_args.suffix))

wcskwds = [ "CRPIX", "CRVAL", "CDELT", "CTYPE", "CUNIT"]
neworder = (1, 2, 0)             # New order of axes (from python's point of view)

for file in files:
    hdu, = pyfits.open(file)
    hdu.data = np.transpose(hdu.data.copy(), neworder)

    orighdr = dict(hdu.header)
    for inew, iold in enumerate(neworder):
        for kwd in wcskwds:
            # FITS uses a convention for the axis order that is inverted from python's order
            hdu.header.update(kwd+str(3-inew), orighdr[kwd+str(3-iold)])
    
    
    hdu.writeto("pv-" + file, clobber=True)


