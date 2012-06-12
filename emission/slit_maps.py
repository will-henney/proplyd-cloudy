import pyfits
import glob
import numpy as np
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="""
Calculate fluxes maps and total emission for a slit from PPV cubes
    """)

parser.add_argument(
    "--suffix", type=str, default="nphi400-nvel200-inc60",  
    help='Last part of cube file names')

parser.add_argument(
    "--onlylinesin", type=str, default=None,
    help='Only process cubes for lines listed in this file (one per line)')

parser.add_argument(
    "--slit", type=float, nargs=4, default=(-1.0, -2.0, 3.0, 2.0),
    metavar=("Xi", "Yi", "Xf", "Yf"),
    help='If present, bounding box of view slit (in units of r0) for the position-position-velocity cube')



cmd_args = parser.parse_args()

cubefiles_found = glob.glob("cube-*-{}.fits".format(cmd_args.suffix))

emlines_found = [ s.split("-")[1].split(".")[0] for s in cubefiles_found ]


if cmd_args.onlylinesin:
    # Read list of selected emission lines if asked for
    with open(cmd_args.onlylinesin) as f:
        selectedlines = f.read().split('\n')
    # But make sure that H beta is in the list
    if not "H__1__4861A" in selectedlines:
        selectedlines.append("H__1__4861A")
    print "Selected lines: ", selectedlines
    # Only calculate lines that are not in the selected list
    emlines = list()
    cubefiles = list()
    for emline, cubefile in zip(emlines_found, cubefiles_found):
        if emline in selectedlines:
            emlines.append(emline)
            cubefiles.append(cubefile)
    print "Lines to be used: ", emlines
else:
    # Otherwise, use all the lines whose emissivity files can be found
    emlines = emlines_found
    cubefiles = cubefiles_found

print "H beta line flux: "
print "Line fluxes relative to H beta = 100:"

for line, cube in zip(emlines, cubefiles):
    hdu, = pyfits.open(cube)

    # construct x-axis slit array
    j0, x0, dx, nx = [hdu.header[kwd] for kwd in "CRPIX1", "CRVAL1", "CDELT1", "NAXIS1"]
    x = x0 + dx*(np.arange(nx) - (j0 - 1))
    xslit = np.where(( x >= Xi ) & ( Xf >= x))[0]
    x = x[:,xmask]
    x = x.reshape((-1, 1, 1))   # make 3D to allow broadcasting with hdu.data

    # construct y-axis slit array
    k0, y0, dy, ny = [hdu.header[kwd] for kwd in "CRPIX2", "CRVAL2", "CDELT2", "NAXIS2"]
    y = y0 + dy*(np.arange(ny) - (k0 - 1))
    yslit = np.where(( y >= Yi ) & ( Yf >= y))[0]
    y = y[:,xmask]
    y = y.reshape((-1, 1, 1))   # make 3D to allow broadcasting with hdu.data


    # Calculate the emission map
    Flux_map = np.sum(hdu.data, axis=0)

    # Calculate the total flux
    Subtotal = np.sum(hdu.data, axis=1)
    Total = np.sum(hdu.data, axis=2)

    hbeta = Total['H__1__4861A']
    print '|'.join(['', line, "%g" % (100.0*Total/hbeta), '']) 

    # Write to 2D FITS files
    for data, id_ in [
        [Flux, "Flux-slit-PPmap"] 
        ]:
        # new HDU to hold this moment map
        newhdu = pyfits.PrimaryHDU(data, hdu.header)
        # remove the kwds associated with the velocity axis
        for kwd in newhdu.header.keys():
            if kwd.endswith("3"):
                del newhdu.header[kwd]
        # save each moment to a file
        newhdu.writeto(id_ + file[4:], clobber=True)


