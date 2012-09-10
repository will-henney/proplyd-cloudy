import pyfits
import glob
import numpy as np
import argparse
import os

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

parser.add_argument(
    "--rotation", type=float, default=0.0,
    help='Angle to the proplyd axis in degrees in the images or spectrum')

cmd_args = parser.parse_args()
cubefiles_found = glob.glob("cube-*-{}.fits".format(cmd_args.suffix))
emlines_found = [ s.split("-")[1].split(".")[0] for s in cubefiles_found ]

xc, yc = cmd_args.center
w, h =  cmd_args.size
xi, yi, xf, yf = xc - 0.5*w, yc - 0.5*h, xc + 0.5*w, yc+0.5*h
alpha = cmd_args.rotation


# Create the folder, if not exist, where the profiles will be
fitsdir = os.path.join(os.getcwd(), 
                       'aperture-xc{:02}-yc{:02}-w{:02}-h{:02}-{}'.format(
        int(10*xc), int(10*yc), int(10*w), int(10*h), cmd_args.suffix))

if not os.path.isdir(fitsdir): 
    os.makedirs(fitsdir)

def saveprofile(data, name, hdr):
    """
    Save a map of some variable to a FITS file
    """
    newhdu = pyfits.PrimaryHDU(data, hdr)
    newhdu.writeto(os.path.join(fitsdir, "profile-{}.fits".format(name)), clobber=True)


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

aperture_profile = dict()
total_aperture_flux = dict()

for line, cube in zip(emlines, cubefiles):

    hdu, = pyfits.open(cube)

    k0, u0, du, nu = [hdu.header[kwd] for kwd in "CRPIX3", "CRVAL3", "CDELT3", "NAXIS3"]

    # construct x-axis slit index
    i0, x0, dx, nx = [hdu.header[kwd] for kwd in "CRPIX1", "CRVAL1", "CDELT1", "NAXIS1"]
    # find indices that correspond to the limits xi and xf
    i1 = i0 + int((xi - x0)/dx)
    i2 = i0 + int((xf - x0)/dx)
    # ensure that we do not exceed the array bounds
    i1 = max(min(i1, nx), 0)
    i2 = max(min(i2, nx), 0)

    # construct y-axis slit index
    j0, y0, dy, ny = [hdu.header[kwd] for kwd in "CRPIX2", "CRVAL2", "CDELT2", "NAXIS2"]
    j1 = j0 + int((yi - y0)/dy)
    j2 = j0 + int((yf - y0)/dy)
    j1 = max(min(j1, ny), 0)
    j2 = max(min(j2, ny), 0)

    cubo = hdu.data

    # rotate the cubes to match the image/spectrum direction

    calpha = np.cos(alpha)
    salpha = np.sin(alpha)

    

    # select the data cube restricted to the slit 
    aperture_cube = cubo[:,j1:j2,i1:i2]

    aperture_profile[line] = np.zeros(nu)
    total_aperture_flux[line] = 0.0

    # add add in to a profile
    aperture_profile[line] = np.sum(np.sum(aperture_cube, axis=-1), axis=-1)

    hdr = pyfits.Header()       # make a new header for the spectrum
    for prefix in ["CRPIX", "CRVAL", "CDELT", "CUNIT", "CTYPE"]:
        # and set the kwds for the first axis from the third axis of the cube
        hdr.update(prefix + "1", hdu.header[prefix + "3"])
    # write the spectrum to a 1D fits image file
    saveprofile(aperture_profile[line], line, hdr)

    # add add add in to a total flux
    total_aperture_flux[line] = aperture_cube.sum()


hbeta = total_aperture_flux["H__1__4861A"]
print "H beta line flux: "
print "Line fluxes relative to H beta = 100:"

for line in emlines:
    flux = total_aperture_flux[line]
    em = line.replace('_', ' ')
    print '|'.join(['', em[:-5], em[-5:], "%g" % (100.0*flux/hbeta), ''])
    
    

    
 



