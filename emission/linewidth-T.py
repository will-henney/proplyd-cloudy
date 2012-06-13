import pyfits
import glob
import numpy as np
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="""
Empirical determination of gas temperature from line moments
    """)


parser.add_argument(
    "--suffix", type=str, default="nphi400-nvel200-inc60",  
    help='Last part of cube file names')

cmd_args = parser.parse_args()

# Options to use either the nebular or auroral linesets
linesets = [
    dict(name="nebular", ha="H__1__6563A", oiii="O__3__5007A", nii="N__2__6584A"),
    dict(name="auroral", ha="H__1__6563A", oiii="TOTL__4363A", nii="N__2__5755A")
    ]

nii_ha_emiss_coeff_ratio = dict(
    nebular=0.5,
    auroral=0.1
    )
dT = 1.e3                       # T difference between N+ and O++ zones


def savemap(data, name):
    """
    Save a map of some variable to a FITS file
    """
    newhdu = pyfits.PrimaryHDU(data)
    newhdu.writeto("LW-{}-{}-{}.fits".format(name, lineset["name"], cmd_args.suffix), clobber=True)


lines = ["ha", "oiii", "nii"]
for lineset in linesets:
    # surface brightness
    S = dict(
        [(line, pyfits.open("M0-sum-{}-{}.fits".format(lineset[line], cmd_args.suffix))[0].data) 
          for line in lines]
        )
    # mean velocity
    V = dict(
        [(line, pyfits.open("M1-mean-{}-{}.fits".format(lineset[line], cmd_args.suffix))[0].data) 
          for line in lines]
        )
    # RMS width
    sig = dict(
        [(line, pyfits.open("M2-sigma-{}-{}.fits".format(lineset[line], cmd_args.suffix))[0].data) 
          for line in lines]
        )

    # First calculate naive T
    T0 = 1.e4*(sig["ha"]**2 - sig["oiii"]**2)/77.34
    savemap(T0, "Naive-T")

    # Second, calculate the correction factors
    
    # Fraction of O++ - Method 1: surface brightness
    f1 = 1.0 - (S["nii"]/S["ha"]) / nii_ha_emiss_coeff_ratio[lineset["name"]] 
    savemap(f1, "f1")

    # Fraction of O++ - Method 2: mean velocities
    f2 = 1.0 - (V["ha"] - V["oiii"])/(V["nii"] - V["oiii"])
    savemap(f2, "f2")

    # Correction factors with the two different ways of estimating f
    for f, fname in (f1, "f1"), (f2, "f2"):
        # Equation (8) of GDHD08
        eta = 1./(1.0 - 0.00957*(1.0 - f))
        # Equation (9) of GDHD08
        epsilon = 0.00957*eta*f*(1.0 - f)*dT/1.e4
        # Equation (10) of GDHD08
        zeta = 0.01293*eta*(1.0 - f)*(sig["nii"]**2 - sig["oiii"]**2)
        # Equation (11) of GDHD08
        xi = 0.01293*eta*f*(1.0 - f)*(V["nii"] - V["oiii"])**2

        savemap(eta, "eta-{}".format(fname))
        savemap(epsilon, "epsilon-{}".format(fname))
        savemap(zeta, "zeta-{}".format(fname))
        savemap(xi, "xi-{}".format(fname))

        # Finally, calculate the corrected temperature
        T = eta*T0 + 1.e4*(epsilon - zeta - xi)
        savemap(T, "Corrected-T-{}".format(fname))



