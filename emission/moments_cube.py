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

cmd_args = parser.parse_args()

files = glob.glob("cube-*-{}.fits".format(cmd_args.suffix))


for file in files:
    hdu, = pyfits.open(file)
    # construct velocity array
    i0, u0, du, n = [hdu.header[kwd] for kwd in "CRPIX3", "CRVAL3", "CDELT3", "NAXIS3"]
    u = u0 + du*(np.arange(n) - (i0 - 1))
    u = u.reshape((-1, 1, 1))   # make 3D to allow broadcasting with hdu.data
    
    # Calculate the velocity moments
    # 0 Moment: Sum (dimensions of hdu.data)
    Sum = np.sum(hdu.data, axis=0)
    # 1 Moment: Mean (dimensions of u)
    Mean = np.sum(hdu.data*u, axis=0) / Sum
    # Remaining moments are central moments about mean
    # 2 Moment: Sigma (dimensions of u)
    Sigma = np.sqrt( np.sum(hdu.data*(u-Mean)**2, axis=0) / Sum )
    # 3 Moment: Skewness (non-dimensional)
    Skew = np.sum(hdu.data*(u-Mean)**3, axis=0) / (Sum * Sigma**3)
    # 4 Moment: Kurtosis (non-dimensional) - subtract 3 so that Gaussian has Kurt = 0
    Kurt = (np.sum(hdu.data*(u-Mean)**4, axis=0) / (Sum * Sigma**4)) - 3.0
    
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
        # save each moment to a file
        newhdu.writeto(id_ + file[4:], clobber=True)
        
        


