"""
Run a series of cloudy proplyd models, automatically adjusting the density


Original version based on ~/Work/Bowshocks/LLobjects/Cloudy/ll1-thon-sd.py

"""

import numpy as np
import scipy
from scipy.optimize import fsolve, bisect
import subprocess
import sys, os
import multiprocessing

from cloudycontrol.model import Model
from cloudycontrol import incident, save, misc, physical

sys.path.append("../src"); import claudia

scriptname = sys.argv[0]

DEGREES = np.pi/180.0           # conversion from degrees to radians

# Proplyd parameters for LV2

def calculate_phi(r0, n0, omega=0.121049705173, alpha=2.6e-13):
    """
    Ionizing flux required for a nominal n0 and r0: Phi = omega r0 n0**2 alpha
    """
    return omega*alpha*r0*n0**2

def calculate_nominal_density(phi, r0, omega=0.121049705173, alpha=2.6e-13, dustfrac=0.4):
    """
    Complementary calculation of density as a function of ionizing flux

    dustfrac is the fraction of EUV photons that are assumed to be absorbed by dust in the ionized gas
    """
    return np.sqrt((1.0 - dustfrac)*phi / (omega*alpha*r0))

def velocity_law_U_of_R(R):
    "Explicit velocity law U(R) requires numerically solving implicit equation"
    def f(U):
	"Function to be zeroed: f(U) = 0"
        return R - velocity_law_R_of_U(U)
    U1, U2 = 1.0, 10.0
    return bisect(f, U1, U2)

def velocity_law_R_of_U(U):
    "Implicit velocity law R(U) from isothermal Dyson model"
    return np.exp(0.25*(U**2-1.0))/np.sqrt(U) 
    
def calculate_hden_Rmax(Rmax, n0):
    """
    Find the density at the illuminated face of the Cloudy model

    This depends on the nominal i-front density, n0, and the radius
    corresponding to z=0 in units of r0, Rmax.

    It also depends on velocity law: n(Rmax) = n0 / Rmax**2 / U(Rmax)

    """
    return n0 /  Rmax**2 / velocity_law_U_of_R(Rmax)
    

def find_max_sound(ccmodel):
    """
    Find the radius at which sound speed is maximum

    Actually we use the temperature to avoid problems with the He ionization front
    """
    # restrict to just the overview file for efficiency
    claudia.SAVETYPES = ["overview"]
    # make a claudia model that corresponds to the cloudycontrol model
    m = claudia.CloudyModel(ccmodel.prefix, indir=ccmodel.indir, outdir=ccmodel.outdir)
    imax_sound = m.pre.cadwind_kms.argmax()
    imax_T = m.ovr.Te.argmax()
    R = args.Rmax - m.ovr.depth/args.r0
    print "Max sound at R=%.4f. Max Te at R=%.4f" % (R[imax_sound], R[imax_T])
    return R[imax_T]



def find_Rsonic_for_hden(hden_Rmax, args, anglepars):
    """
    Runs a Cloudy model to find the dimensionless radius R of the sonic point
    
    Principal argument: hden_Rmax

    Other vars are taken from args (command line parameters) and anglepars

    it - problematic when called as part of root finding routine
    """
    cloudy_cmd = ["time", args.cloudyexec] 
    # files are simply named after iteration number and density
    filename = "it%.2in%.3e" % (it, hden_Rmax)
    print "Calculating " + os.path.join(anglepars.datapath, filename)
    m = Model(filename, dryrun=False, indir=anglepars.datapath, outdir=anglepars.datapath, 
              verbose=True, cloudy_cmd=cloudy_cmd)
    if args.fastcloudy:
        dust = "FastOrion"
    else:
        dust = "Orion"
    m.write(physical.proplyd(args.r0, hden_Rmax, Rmax=args.Rmax, 
                             W=args.W, x0=args.x0, dust=dust))
    m.write(misc.optimize())
    m.write(misc.stopping())
    if args.fastcloudy:
        m.write(misc.iterate(iterations=0))
    else:
        m.write(misc.iterate())
    m.write(incident.background())
    m.write(incident.star(anglepars.logPhiH, args.Tstar, atmosphere=args.atmosphere, log_g=4.1)) # th1C
    m.write(save.default())
    if args.savecont:
        m.write("""\
save every last ionizing continuum 1.0 0 ".icont"
""")
    m.run()
    # Assume that sonic point is at sound speed maximum (actually we use Te)
    return find_max_sound(m)





# Automatic adjustment to try and fix i-front
beta_index = -4.5


class AnglePars(object):
    pass


if __name__ == '__main__':
    import warnings
    import argparse

    warnings.filterwarnings("ignore") # we shouldn't really do this


    parser = argparse.ArgumentParser(
        description="Run a series of cloudy proplyd models, automatically adjusting the density",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument("--cloudyexec", "-c", type=str, default='cloudy.exe', 
                        help='Name of Cloudy executable')
    parser.add_argument("--ntheta", "-n", type=int, default=3,
                        help="Number of incident radiation angles to calculate between 0 and 90 deg")
    parser.add_argument("--rootdir", "-r", type=str, default="proplyd-auto-models", 
                        help='Directory where output will be written')
    parser.add_argument("--tolerance", "-t", type=float, default=1.e-3,
                        help="Tolerance to which Rsonic=1 must be satisfied for root finding")
    parser.add_argument("--r0", type=float, default=8e14,
                        help="Proplyd sonic radius (in cm)")
    parser.add_argument("--Rmax", type=float, default=9.0,
                        help="Dimensionless outer radius of Cloudy models")
    parser.add_argument("--W", type=float, default=30.0,
                        help="Thickness of i-front in units of mfp at 1 Rydberg")
    parser.add_argument("--x0", type=float, default=10.8,
                        help="Offset between sonic point and i-front in units of W/2")
    parser.add_argument("--diffuseBeta", type=float, default=0.004,
                        help="Diffuse EUV flux to stellar flux ratio")
    parser.add_argument("--logPhiH", type=float, default=14.54,
                        help="Log10 of stellar EUV flux")
    parser.add_argument("--atmosphere", type=str, default="WM", choices=["WM", "BB", "TL"],
                        help="Type of stellar atmosphere model")
    parser.add_argument("--Tstar", type=float, default=3.7e4,
                        help="Stellar effective temperature in K")
    parser.add_argument("--multiprocess", "-m", type=bool, default=False,
                        help="Whether to use multiprocessing (currently not implemented)")
    parser.add_argument("--fastcloudy", "-f", type=bool, default=False,
                        help="Whether to make cloudy run faster (at the expense of accuracy)")
    parser.add_argument("--savecont", "-s", type=bool, default=False,
                        help="Whether to make cloudy save the continuum for every zone")
    args = parser.parse_args()

    
    
    if args.multiprocess:
        raise NotImplementedError # Delete this when it is implemented properly
        pool = multiprocessing.Pool() # by default makes one process per core

    # ID for this proplyd model
    modelid = "%s%.6i-phi%.2f-r%.2f" % (args.atmosphere, args.Tstar, args.logPhiH, np.log10(args.r0))
    # the ID may have to be refined by adding more parameters later

    thetas = np.linspace(0.0, 90.0, args.ntheta)
    for theta in thetas:
        anglepars = AnglePars() # container for miscellaneous parameters
        thetaid = "th%.2i" % (theta)
        # All the files for a given angle now go in the same directory...
        anglepars.datapath = os.path.join(args.rootdir, modelid, thetaid)
        PhiH = (10**args.logPhiH) * (np.cos(theta * DEGREES) + args.diffuseBeta)
        n0 = calculate_nominal_density(PhiH, args.r0)
        anglepars.logPhiH = np.log10(PhiH)
        hden_Rmax = calculate_hden_Rmax(args.Rmax, n0) # first guess for outer density

        converged = False
        it = 0
        # bestHi and bestLo represent the current brackets on the density
        # Initially we take very generous brackets 
        bestHi = 10.0*hden_Rmax
        bestLo = 0.1*hden_Rmax
        print "Initial brackets: [%.2e, %.2e]" % (bestLo, bestHi) 
        while not converged:
            it += 1
            Rsonic = find_Rsonic_for_hden(hden_Rmax, args, anglepars)
            # Update bounds on solution
            # This assumes that d R / d hden is positive, which should always be the case
            if Rsonic > 1.0:
                bestHi = min(bestHi, hden_Rmax)
            else:
                bestLo = max(bestLo, hden_Rmax)
            print "New brackets: [%.2e, %.2e]" % (bestLo, bestHi) 
            if abs(Rsonic - 1.0) > args.tolerance:
                print "Sonic point is in wrong place (R=%.4f) - adjust density and try again" % (Rsonic)
                # IF NOT, FIND CORRECTED OUTER DENSITY
                if it <= 2:
                    # pivot the first couple of times
                    hden_Rmax *= Rsonic**beta_index
                    print "Pivoting ..."
                else:
                    # then use bisection
                    hden_Rmax = 0.5*(bestHi + bestLo)
                    print "Bisecting ..."

            else:
                converged = True

        print "Angle %i converged in %i iteration%s" % (theta, it, "s" if it > 1 else "")

    if args.multiprocess:
        # Clean up everything
        pool.close()
        pool.join()


