import pyfits
import glob
import numpy as np
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="""
Calculate velocity moments of PPV cubes
    """)

parser.add_argument(
    "--suffix", type=str, default="nphi400-nvel200-inc60",  
    help='Last part of cube file names')

parser.add_argument(
    "--onlylinesin", type=str, default=None,
    help='Only process cubes for lines listed in this file (one per line)')

parser.add_argument(
    "--rotation", type=float, default=00.0,
    help='Angle to the proplyd axis in degrees in the images or spectrum')

parser.add_argument(
    "--pixelsize", type=float, default=0.1, 
    help='Size of pixels (in units of r0) for the moments maps')

parser.add_argument(
    "--center", type=float, nargs=2, default=(1.0, 0.0),
    metavar=("x", "y"),
    help='Coordinates of the center of the aperture')

parser.add_argument(
    "--aperturesize", type=float, nargs=2, default=(1.0, 1.0), metavar=("w", "h"),
    help='Size of aperture (width, height)')

cmd_args = parser.parse_args()

cubefiles_found = glob.glob("cube-*-{}.fits".format(cmd_args.suffix))

emlines_found = [ s.split("-")[1].split(".")[0] for s in cubefiles_found ]

xc, yc = cmd_args.center
w, h =  cmd_args.size

alpha = cmd_args.rotation
ca = np.cos(np.radians(alpha))
sa = np.sin(np.radians(alpha))

pixely = pixelx = cmd_args.pixelsize

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

# Create the folder, if not exist, where the profiles will be
fitsdir = os.path.join(os.getcwd(), 
                       'moments-xc{:02}-yc{:02}-w{:02}-h{:02}-rot{:02}-{}'.format(
        int(10*xc), int(10*yc), int(10*w), int(10*h), int(alpha), cmd_args.suffix))

# Construct the aperture
    

for file in cubefiles:
    hdu, = pyfits.open(file)
    # construct velocity array
    i0, u0, du, n = [hdu.header[kwd] for kwd in "CRPIX3", "CRVAL3", "CDELT3", "NAXIS3"]
    u = u0 + du*(np.arange(n) - (i0 - 1))
    u = u.reshape((-1, 1, 1))   # make 3D to allow broadcasting with hdu.data
    # construct x-axis slit index
    i0, x0, dx, nx = [hdu.header[kwd] for kwd in "CRPIX1", "CRVAL1", "CDELT1", "NAXIS1"]
    # construct y-axis slit index
    j0, y0, dy, ny = [hdu.header[kwd] for kwd in "CRPIX2", "CRVAL2", "CDELT2", "NAXIS2"]

    # construct an x array
    xf = x0 + (nx*dx)
    x = np.linspace(x0,xf,nx)

    # construct an y array
    yf = y0 + (ny*dy)
    y = np.linspace(y0,yf,ny)

    # Promote coordinate axes to 2d arrays
    X, Y = np.meshgrid(x, y)

    # Rotate to the frame of the aperture
    P = (X-xc)*ca - (Y-yc)*sa
    Q = (X-xc)*sa + (Y-yc)*ca

    mask = (abs(P) < 0.5*w) & (abs(Q) < 0.5*h)

    # select the data cube restricted to the slit 
    aperture_cube = np.where( mask, hdu.data, 0.0 )

    # Calculate the velocity moments
    # 0 Moment: Sum (dimensions of hdu.data)
    Sum = np.sum(aperture_cube, axis=0)
    # 1 Moment: Mean (dimensions of u)
    Mean = np.sum(aperture_cube*u, axis=0) / Sum
    # Remaining moments are central moments about mean
    # 2 Moment: Sigma (dimensions of u)
    Sigma = np.sqrt( np.sum(aperture_cube*(u-Mean)**2, axis=0) / Sum )
    # 3 Moment: Skewness (non-dimensional)
    Skew = np.sum(aperture_cube*(u-Mean)**3, axis=0) / (Sum * Sigma**3)
    # 4 Moment: Kurtosis (non-dimensional) - subtract 3 so that Gaussian has Kurt = 0
    Kurt = (np.sum(aperture_cube*(u-Mean)**4, axis=0) / (Sum * Sigma**4)) - 3.0
    
    # Write to 2D FITS files
    for data, id_ in [
        [Sum, "M0-sum"], 
        [Mean, "M1-mean"], 
        [Sigma, "M2-sigma"], 
        [Skew, "M3-skewness"], 
        [Kurt, "M4-kurtosis"] 
        ]:
        # new HDU to hold this moment map
        newhdu = pyfits.PrimaryHDU(data, hdu.header)
        # remove the kwds associated with the velocity axis
        for kwd in newhdu.header.keys():
            if kwd.endswith("3"):
                del newhdu.header[kwd]
        # save each moment to a file in the folder        
        newhdu.writeto(os.path.join(fitsdir, id_ + file[4:], clobber=True)


