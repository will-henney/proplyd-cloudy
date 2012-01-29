"""
RGB plots of emissivity in the Ne-Te plane

Color represents different emission lines

"""

import bivar                            # my bivariate distribution plotting library
import numpy as np
import os, sys, argparse, glob

sys.path.append("../../../src")       # make sure we can find claudia.py
import claudia

# Avoid verbose error messages from numpy during the reading of the Cloudy files
import warnings, argparse
warnings.filterwarnings("ignore")

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
    )
parser.add_argument(
    "--lineset", type=str, default="SNO", choices=linesets.keys(),
    help='Which set of 3 emission lines to use as R, G, B channels')

cmd_args = parser.parse_args()

emlines = [s.replace(' ', '_') for s in linesets[cmd_args.lineset]]

thetadirs =  [
    os.path.basename(p) for p in 
    glob.glob(os.path.join(cmd_args.modeldir, "th??"))
    ]

# Lists for accumulating the data
Eden_list = []
Te_list = []
Red_weight_list = []
Green_weight_list = []
Blue_weight_list = []

Theta = np.array([np.radians(float(s[2:])) for s in thetadirs])
NJ = len(Theta)
for j, thetadir in enumerate(thetadirs):
    # Read in the Cloudy model for each angle
    claudia.CloudyModel.indir = os.path.join(cmd_args.modeldir, thetadir)
    claudia.CloudyModel.outdir = claudia.CloudyModel.indir
    # find the last iteration
    modelinputpath = glob.glob(os.path.join(claudia.CloudyModel.indir, "*.in"))[-1]
    # extract only the basename of the file, without .in extension
    modelname = os.path.splitext(os.path.basename(modelinputpath))[0]
    print "Indir: ", claudia.CloudyModel.indir
    print "Modelname: ", modelname
    m = claudia.CloudyModel(modelname) 
    
    # Calculate required variables as function of radius for this theta
    R = cmd_args.Rmax - m.ovr.depth/cmd_args.r0 # array of dimensionless radius
    Eden = 10**(m.ovr.eden)
    Te = 10**(m.ovr.Te)
    Emissivities = [10**(m.em[emline]) for emline in emlines]

    ctheta = np.cos(Theta[j])
    stheta = np.sin(Theta[j])
    jneg = max(0, j - 1)
    jpos = min(NJ-1, j + 1)
    dmu = -0.5*(np.cos(Theta[jpos]) - np.cos(Theta[jneg]))

    # Use a simple accumulator loop, even though it is inefficient
    NI = len(R)
    for i in range(NI):
        ineg = max(0, i - 1)
        ipos = min(NI-1, i + 1)
        dr = 0.5*(R[ipos] - R[ineg])
        dvol =  2.* np.pi * dmu * (R[i]**2) * dr
        Eden_list.append(Eden[i])
        Te_list.append(Te[i])
        Red_weight_list.append(Emissivities[0][i]*dvol)
        Green_weight_list.append(Emissivities[1][i]*dvol)
        Blue_weight_list.append(Emissivities[2][i]*dvol)



Eden_arr = np.array(Eden_list).astype(float)
Te_arr = np.array(Te_list).astype(float)
weights = [
    np.array(Red_weight_list).astype(float),
    np.array(Green_weight_list).astype(float),
    np.array(Blue_weight_list).astype(float)
    ]

# Now make the plots
import pyx
pyx.text.set(mode="latex")
# pyx.text.preamble("""\usepackage{mathptmx}""")
bivar.printextra = 0
bivar.Graph.figwidth = 6
bivar.Graph.figheight = 6
bivar.PlotVariable.n = 20                               # size of pdf images

Eden_var = bivar.PlotVariable(np.log10(Eden_arr))
Eden_var.setminmaxn(min=3.0, max=6.0)
Eden_var.settitle(r'Electron density, \(\log_{10} (N_\mathrm{e}/\mathrm{cm}^{-3})\)', 'Eden')

Te_var = bivar.PlotVariable(Te_arr)
Te_var.setminmaxn(min=5000.0, max=12000.)
Te_var.settitle(r'Gas temperature, \(T\), K', 'Te')

g = bivar.Graph(Eden_var, Te_var, weights=weights, gamma=3.0, statslevel=1)


g.writePDFfile('NT-plane-%s' % (cmd_args.lineset))



