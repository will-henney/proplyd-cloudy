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
    metavar=("xi", "yi", "xf", "yf"),
    help='If present, bounding box of view slit (in units of r0) for the position-position-velocity cube')


cmd_args = parser.parse_args()
cubefiles_found = glob.glob("cube-*-{}.fits".format(cmd_args.suffix))
emlines_found = [ s.split("-")[1].split(".")[0] for s in cubefiles_found ]


xi, yi, xf, yf = cmd_args.slit


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

yv_slit_map = dict()
slit_profile = dict()
total_slit_flux = dict()

for line, cube in zip(emlines, cubefiles):

    hdu, = pyfits.open(cube)

    i0, u0, du, nu = [hdu.header[kwd] for kwd in "CRPIX3", "CRVAL3", "CDELT3", "NAXIS3"]

    # construct x-axis slit index
    j0, x0, dx, nx = [hdu.header[kwd] for kwd in "CRPIX1", "CRVAL1", "CDELT1", "NAXIS1"]
    x = x0 + dx*(np.arange(nx) - (j0 - 1))
    xslit = np.where(( x >= xi ) & ( xf >= x))[0]
    XI = min(xslit)
    XF = max(xslit)
    NX = len(xslit)

    # construct y-axis slit index
    k0, y0, dy, ny = [hdu.header[kwd] for kwd in "CRPIX2", "CRVAL2", "CDELT2", "NAXIS2"]
    y = y0 + dy*(np.arange(ny) - (k0 - 1))
    yslit = np.where(( y >= yi ) & ( yf >= y))[0]
    YI = min(yslit)
    YF = max(yslit)
    NY = len(yslit)

    cubo = hdu.data
    # select the data cube restricted to the slit 
    cubito = cubo[:,XI:XF,YI:YF]

    
    yv_slit_map[line] = np.zeros((nu,NY))
    slit_profile[line] = np.zeros(nu)
    total_slit_flux[line] = 0.0


    # add in to a pv map
    yv_slit_map[line] = np.sum(cubito, axis=1)
    # add add in to a profile
    slit_profile[line] = np.sum(yv_slit_map[line], axis=-1)
    # add add add in to a total flux
    total_slit_flux[line] = cubito.sum()

    print yv_slit_map[line].shape, slit_profile[line].shape, total_slit_flux[line].shape

hbeta = total_slit_flux["H__1__4861A"]
print "H beta line flux: "
print "Line fluxes relative to H beta = 100:"

for line in emlines:
    flux = total_slit_flux[line]
    em = line.replace('_', ' ')
    print '|'.join(['', em[:-5], em[-5:], "%g" % (100.0*flux/hbeta), ''])
    
    

    
 



