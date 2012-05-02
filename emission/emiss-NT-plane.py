"""
RGB plots of emissivity in the Ne-Te plane

Color represents different emission lines

"""

import bivar                            # my bivariate distribution plotting library
import numpy as np
import os, sys, argparse, glob
import pyfits

# Parse command line arguments
parser = argparse.ArgumentParser(
    description="RGB emissivity in Ne-Te plane for a Cloudy proplyd model")

parser.add_argument(
    "--modeldir", type=str, default=".",  
    help='Path to proplyd directory model')
parser.add_argument(
    "--r0", type=float,
    help='Radius of proplyd ionization front in cm (no default so we do not forget)')
parser.add_argument(
    "--Rmax", type=float, default=9.0,
    help='Radius of outer face of cloudy model in units of r0')

linesets = dict(
    SNO = ["S II  6731A", "N  2  6584A", "O  3  5007A"],
    SNOa = ["S II  4070A", "N  2  5755A", "TOTL  4363A"],
    ONO = ["O  1  6300A", "N  2  6584A", "O  3  5007A"],
    NHO = ["N  2  6584A", "H  1  6563A", "O  3  5007A"],
    SSS = ["S II  6716A", "S II  6731A", "S II  4070A"],
    SONe = ["S  3  6312A", "O  3  5007A", "Ne 3  3869A"],
    )
parser.add_argument(
    "--lineset", type=str, default="SNO", choices=linesets.keys(),
    help='Which set of 3 emission lines to use as R, G, B channels')
parser.add_argument(
    "--rebin", type=str, default="linear-101",
    help="Which set of rebinned data to use")
parser.add_argument(
    "--N", type=int, default=50,
    help="Number of pixels (NxN) in Ne-Te image")
parser.add_argument(
    "--gamma", type=float, default=2.0,
    help="Gamma correction for Ne-Te image")
parser.add_argument(
    "--mode", type=str, default='RGB', choices=['RGB', 'felt tips', 'oil paints'],
    help="How to combine the three color channels")
parser.add_argument(
    "--redpen", type=str, default='#f06', 
    help="Color for red pen in 'felt tips' and 'oil paints' modes")
parser.add_argument(
    "--greenpen", type=str, default='#5f5', 
    help="Color for green pen in 'felt tips' and 'oil paints' modes")
parser.add_argument(
    "--bluepen", type=str, default='#18f', 
    help="Color for blue pen in 'felt tips' and 'oil paints' modes")



cmd_args = parser.parse_args()

emlines = [s.replace(' ', '_') for s in linesets[cmd_args.lineset]]


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

# Read in interpolated data from the FITS files
Eden = 10**(read_fits("ovr-eden"))
Te = 10**(read_fits("ovr-Te"))
Emissivities = [10**read_fits("em-%s" % (emline)) for emline in emlines]

R = read_fits("R")
dR = np.diff(R)
# Repeat last column of dR so that it has same shape as other arrays
dR = np.hstack((dR, dR[:,-1:]))

mu = read_fits("mu")
# Find number of mu points by parsing name of rebin subdir
NMU = int(cmd_args.rebin.split('-')[-1])
dmu = 1.0/(NMU-1)
dvol = 2.0*np.pi * dmu * R**2 * dR

weights = [em*dvol for em in Emissivities]

# Now make the plots
import pyx
pyx.text.set(mode="latex")
# pyx.text.preamble("""\usepackage{mathptmx}""")
bivar.printextra = 0
bivar.Graph.figwidth = 6
bivar.Graph.figheight = 6
bivar.PlotVariable.n = cmd_args.N # size of pdf images

Eden_var = bivar.PlotVariable(np.log10(Eden))
Eden_var.setminmaxn(min=3.0, max=6.0)
Eden_var.settitle(r'Electron density, \(\log_{10} (n_\mathrm{e}/\mathrm{cm}^{-3})\)', 'Eden')

Te_var = bivar.PlotVariable(Te)
Te_var.setminmaxn(min=5000.0, max=12000.)
Te_var.settitle(r'Gas temperature, \(T\), K', 'Te')

g = bivar.Graph(Eden_var, Te_var, weights=weights, gamma=cmd_args.gamma, 
                statslevel=1, composition_mode=cmd_args.mode, 
                pen_colors = [cmd_args.redpen, cmd_args.greenpen, cmd_args.bluepen],
                channel_order = [2, 0, 1] # blue, red, green
                # channel_order = [2, 1, 0] # blue, green, red
                )



g.writePDFfile('NT-plane-%s' % (cmd_args.lineset))



