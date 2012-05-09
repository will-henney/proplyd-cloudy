"""
Take the models and rebin them at a finer grid in theta

Write out a FITS file for each variable of each Cloudy output file
"""
import numpy as np
import os, sys, argparse, glob
import scipy
import scipy.interpolate
import pyfits

sys.path.append("../../../src")       # make sure we can find claudia.py
import claudia

# Avoid verbose error messages from numpy during the reading of the Cloudy files
import warnings, argparse
warnings.filterwarnings("ignore")

# Parse command line arguments
parser = argparse.ArgumentParser(
    description="Interpolate proplyd models to finer angle grid (and uniform z grid)")

parser.add_argument(
    "--modeldir", type=str, default=".",  
    help='Path to proplyd directory model')
parser.add_argument(
    "--r0", type=float,
    help='Radius of proplyd ionization front in cm (no default so we do not forget)')
parser.add_argument(
    "--Rmax", type=float, default=9.0,
    help='Radius of outer face of cloudy model in units of r0')
parser.add_argument(
    "--nangles", type=int, default=101,
    help='Number of interpolation angles to use (interpolation is actually in mu = cos theta)')
parser.add_argument(
    "--method", type=str, default='linear', choices=['nearest', 'linear', 'cubic'],
    help='The interpolation method could be linear, cubic or nearest neighbour')

cmd_args = parser.parse_args()

def get_cloudy_model(thetadir):
    claudia.CloudyModel.indir = os.path.join(cmd_args.modeldir, thetadir)
    claudia.CloudyModel.outdir = claudia.CloudyModel.indir
    # find the last iteration
    modelinputpath = glob.glob(os.path.join(claudia.CloudyModel.indir, "*.in"))[-1]
    # extract only the basename of the file, without .in extension
    modelname = os.path.splitext(os.path.basename(modelinputpath))[0]
    print "Indir: ", claudia.CloudyModel.indir
    print "Modelname: ", modelname
    return claudia.CloudyModel(modelname) 
    

# First read all the models into memory
thetadirs =  [
    os.path.basename(p) for p in 
    glob.glob(os.path.join(cmd_args.modeldir, "th??"))
    ]

# values of mu of the individual models

mu_pts = []                     # list of mus for all depths of all models
R_pts = []                      # list of depths for all depths of all models


# Use the first angle (th=0) as template for setting up the depth
# array and the data structures
template_model = get_cloudy_model(thetadirs[0])

##
## Set up the data structure for the Cloudy dependent variables
##

# Which cloudy save files to grab the data columns from
tablenames = ["ovr", "em", "ion_c", "ion_n", "ion_o", "ion_s", "phy", "pdr"]

# datatables is a dictionary of tables, with each item corresponding
# to a Cloudy save file (e.g., ovr, em, pre, etc)
datatables = dict()             
for tabkey in tablenames: 
    # Each table item is a dictionary of variables (columns in the save file)
    datatables[tabkey] = dict()
    columnkeys = template_model.__dict__[tabkey].dtype.names[1:] # skip the first column since it is depth
    for columnkey in columnkeys:
        # Each variable starts as an empty list, which will be filled
        # in later with the values from all the angles
        datatables[tabkey][columnkey] = []

##
## For all the model angles, fill in the independent variable lists
## (R_pts, mu_pts) and all the dependent variable lists (datavalues)
##
for j, thetadir in enumerate(thetadirs):
    # Read in the Cloudy model for each angle
    m = get_cloudy_model(thetadir)
    # convert theta to float then find mu = cos(th)
    mu = np.cos(np.radians(float(thetadir[2:]))) 
    # number of grid points in this cloudy model
    N = len(m.ovr.depth)
    R = cmd_args.Rmax - m.ovr.depth/cmd_args.r0
    R_pts.extend(R) # stuff the radii into the long list of points
    mu_pts.extend([mu]*N) # repeat the mu value N times
    # Now, for every model put each column onto end of the correct datavalues list
    for tabkey, table in datatables.items():
        for columnkey, datavalues in table.items():
            datavalues.extend(m.__dict__[tabkey][columnkey])
  
# Set up array of points (R, mu) where we know the data values
points = np.array(zip(R_pts, mu_pts))
print points

# We want to interpolate onto a grid of mu (uniform) and R (non-uniform)

# # For the common depth scale, we just take that of the first model (th=0)
# R_grid = cmd_args.Rmax - template_model.ovr.depth/cmd_args.r0 # dimensionless radius
# Change to use all the points that ever were
R_grid = np.array(sorted(R_pts))
# Angle grid is uniform in cos(th)
mu_grid = np.linspace(0.0, 1.0, cmd_args.nangles)

# Make them both 2-dimensional
R_grid, mu_grid = np.meshgrid(R_grid, mu_grid)

# This is where we will put all the output FITS files
fitsdir = os.path.join(cmd_args.modeldir, 
                       'rebin-%s-%i' % (cmd_args.method, cmd_args.nangles))
if not os.path.isdir(fitsdir): 
    os.makedirs(fitsdir)

def write_fits(data, fileid):
    """Save array DATA to file FITSDIR/FILEID.fits"""
    hdu = pyfits.PrimaryHDU(data)
    hdu.writeto(os.path.join(fitsdir, fileid + '.fits'), clobber=True)

# Save each independent variable grid to a file
write_fits(R_grid, 'R')
write_fits(mu_grid, 'mu')

# Carry out the interpolation for each dependent variable and save a file
for tabkey, table in datatables.items(): 
    for columnkey, datavalues in table.items():
        datagrid = scipy.interpolate.griddata(
            points, np.array(datavalues), 
            (R_grid, mu_grid), 
            method=cmd_args.method
            )
        write_fits(datagrid, "%s-%s" % (tabkey, columnkey))
        



