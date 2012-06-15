import pyfits
import glob
import numpy as np
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="""
Calculate line intensities through an aperture (all dimensions are in units of r0)
    """)

parser.add_argument(
    "--suffix", type=str, default="nphi400-nvel200-inc60",  
    help='Last part of cube file names')

parser.add_argument(
    "--onlylinesin", type=str, default=None,
    help='Give intensities for only those lines listed in this file (one per line)')

parser.add_argument(
    "--center", type=float, nargs=2, default=(1.0, 0.0),
    metavar=("x", "y"),
    help='Coordinates of the center of the aperture')

parser.add_argument(
    "--size", type=float, nargs=2, default=(1.0, 1.0), metavar=("w", "h"),
    help='Size of aperture (width, height)')


cmd_args = parser.parse_args()
cubefiles_found = glob.glob("cube-*-{}.fits".format(cmd_args.suffix))
emlines_found = [ s.split("-")[1].split(".")[0] for s in cubefiles_found ]


xc, yc = cmd_args.center
w, h =  cmd_args.size
xi, yi, xf, yf = xc - 0.5*w, yc - 0.5*h, xc + 0.5*w, yc+0.5*h

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

    k0, u0, du, nu = [hdu.header[kwd] for kwd in "CRPIX3", "CRVAL3", "CDELT3", "NAXIS3"]

    # construct x-axis slit index
    i0, x0, dx, nx = [hdu.header[kwd] for kwd in "CRPIX1", "CRVAL1", "CDELT1", "NAXIS1"]
    # find indices that correspond to the limits xi and xf
    i1 = i0 + int((xi - x0)/dx)
    i2 = i0 + int((xf - x0)/dx)
    # ensure that we do neot exceed the array bounds
    i1 = max(min(i1, nx), 0)
    i2 = max(min(i2, nx), 0)


    # construct y-axis slit index
    j0, y0, dy, ny = [hdu.header[kwd] for kwd in "CRPIX2", "CRVAL2", "CDELT2", "NAXIS2"]
    j1 = j0 + int((yi - y0)/dy)
    j2 = j0 + int((yf - y0)/dy)
    j1 = max(min(j1, ny), 0)
    j2 = max(min(j2, ny), 0)

    cubo = hdu.data
    # select the data cube restricted to the slit 
    cubito = cubo[:,j1:j2,i1:i2]

    
    slit_profile[line] = np.zeros(nu)
    total_slit_flux[line] = 0.0


    # add add in to a profile
    slit_profile[line] = np.sum(np.sum(cubito, axis=-1), axis=-1)
    # add add add in to a total flux
    total_slit_flux[line] = cubito.sum()


hbeta = total_slit_flux["H__1__4861A"]
print "H beta line flux: "
print "Line fluxes relative to H beta = 100:"

for line in emlines:
    flux = total_slit_flux[line]
    em = line.replace('_', ' ')
    print '|'.join(['', em[:-5], em[-5:], "%g" % (100.0*flux/hbeta), ''])
    
    

    
 



