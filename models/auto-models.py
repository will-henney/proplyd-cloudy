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
SIGMA0 = 6.30e-18               # Threshold photoionization x-section

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

def find_max_sound(ccmodel, use_T_not_c=False):
    """
    Find the radius at which sound speed is maximum

    Option to use the temperature to avoid problems with the He ionization front
    """
    # restrict to just the overview file for efficiency
    claudia.SAVETYPES = ["overview"]
    # make a claudia model that corresponds to the cloudycontrol model
    m = claudia.CloudyModel(ccmodel.prefix, indir=ccmodel.indir, outdir=ccmodel.outdir)
    imax_sound = m.pre.cadwind_kms.argmax()
    imax_T = m.ovr.Te.argmax()
    R = args.Rmax - m.ovr.depth/args.r0
    print "Max sound at R=%.4f. Max Te at R=%.4f" % (R[imax_sound], R[imax_T])
    if use_T_not_c:
        return R[imax_T]
    else:
        return R[imax_sound]

def find_ifront(ccmodel):
    """
    Find the radius at which H ionization fraction = 0.5
    """
    # restrict to just the overview file for efficiency
    claudia.SAVETYPES = ["overview"] # why does this not do anything?
    # make a claudia model that corresponds to the cloudycontrol model
    m = claudia.CloudyModel(ccmodel.prefix, indir=ccmodel.indir, outdir=ccmodel.outdir)
    ifrac = 10**m.ovr.HII
    R = args.Rmax - m.ovr.depth/args.r0
    Rif = find_x_at_y(0.5, R, ifrac)
    return Rif


def find_x_at_y(ystar, x, y):
    """
    Find the xstar where y=ystar, using a piecewise linear approximation 
    to discrete vectors x and y (of equal lengths) 
    """
    if y[0] > ystar and y[-1] < ystar:
        # case of y falling with index
        i2 = len(y[y>ystar])
    elif y[0] < ystar and y[-1] > ystar:
        # case of y increasing with index
        i2 = len(y[y<ystar])
    else:
        raise ValueError, "y array must include the target value %.3f" % (ystar)
    i1 = i2 - 1
    y1, y2 = y[i1], y[i2]
    x1, x2 = x[i1], x[i2]
    xstar = x1 + (ystar-y1)*(x2-x1)/(y2-y1)
    return xstar

def setup_and_run_cloudy_model(args, extras):
    cloudy_cmd = ["time", args.cloudyexec] 
    # files are simply named after iteration number and density
    filename = "it%.2in%.3ex0%.2f" % (extras.it, extras.hden, extras.x0)
    print "Calculating " + os.path.join(extras.datapath, filename)
    m = Model(filename, dryrun=False, indir=extras.datapath, outdir=extras.datapath, 
              verbose=True, cloudy_cmd=cloudy_cmd)
    if args.fastcloudy:
        dust = "FastOrion"
    else:
        dust = "Orion"
    m.write(physical.proplyd(args.r0, extras.hden, Rmax=args.Rmax, 
                             W=args.W, x0=extras.x0, dust=dust))
    m.write(misc.optimize())
    m.write(misc.stopping())
    if args.fastcloudy:
        m.write(misc.iterate(iterations=0))
    else:
        m.write(misc.iterate())
    m.write(incident.background())
    m.write(incident.star(extras.logPhiH, args.Tstar, atmosphere=args.atmosphere, log_g=4.1)) # th1C
    m.write(save.default())
    if args.savecont:
        m.write("""\
save every last ionizing continuum 1.0 0 ".icont"
""")
    m.run()
    return m
    






def find_Rsonic_for_hden(hden_Rmax, args, extras):
    """
    Runs a Cloudy model to find the dimensionless radius R of the sonic point
    
    Principal argument: hden_Rmax

    Other vars are taken from args (command line parameters) and extras

    it - problematic when called as part of root finding routine
    """
    extras.hden = hden_Rmax
    extras.it += 1
    print "Running cloudy with hden = %.3e" % (hden_Rmax)
    m = setup_and_run_cloudy_model(args, extras)
    # Assume that sonic point is at sound speed maximum (actually we use Te)
    return find_max_sound(m) - 1.0

def find_Rif_for_x0(x0, args, extras):
    """
    Runs a Cloudy model to find the dimensionless radius R of the i-front
    
    Principal argument: x0

    Other vars are taken from args (command line parameters) and extras

    it - problematic when called as part of root finding routine
    """
    extras.x0 = x0
    dR = (0.5*args.W*extras.x0) / \
        (args.r0*SIGMA0*extras.hden*args.Rmax*args.Rmax*velocity_law_U_of_R(args.Rmax))
    extras.it += 1
    print "Running cloudy with x0 = %.2f, dR = %.5f" % (x0, dR)
    m = setup_and_run_cloudy_model(args, extras)
    Rif = find_ifront(m)
    print "I-front radius = %.4f, target = %.4f" % (Rif, (1.0 - dR))
    return Rif - (1.0 - dR)





# Automatic adjustment to try and fix i-front
beta_index = -4.5


class ExtraPars(object):
    pass


if __name__ == '__main__':
    import warnings
    import argparse
    import scipy.optimize

    warnings.filterwarnings("ignore") # we shouldn't really do this


    parser = argparse.ArgumentParser(
        description="Run a series of cloudy proplyd models, automatically adjusting the density",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument("--cloudyexec", type=str, default='cloudy.exe', 
                        help='Name of Cloudy executable')
    parser.add_argument("--ntheta", "-n", type=int, default=3,
                        help="Number of incident radiation angles to calculate between 0 and 90 deg")
    parser.add_argument("--rootdir", type=str, default="proplyd-auto-models", 
                        help='Directory where output will be written')
    parser.add_argument("--tolerance", type=float, default=1.e-3,
                        help="Tolerance to which Rsonic=1 must be satisfied for root finding")
    parser.add_argument("--r0", type=float, default=8e14,
                        help="Proplyd sonic radius (in cm)")
    parser.add_argument("--Rmax", type=float, default=9.0,
                        help="Dimensionless outer radius of Cloudy models")
    parser.add_argument("--W", type=float, default=30.0,
                        help="Thickness of i-front in units of mfp at 1 Rydberg")
    parser.add_argument("--initialx0", type=float, default=10.8,
                        help="Guess at offset between sonic point and i-front in units of W/2 (this parameter is adjusted until the models converge)")
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
        extras = ExtraPars() # container for miscellaneous parameters
        thetaid = "th%.2i" % (theta)
        # All the files for a given angle now go in the same directory...
        extras.datapath = os.path.join(args.rootdir, modelid, thetaid)
        PhiH = (10**args.logPhiH) * (np.cos(theta * DEGREES) + args.diffuseBeta)
        n0 = calculate_nominal_density(PhiH, args.r0)
        extras.logPhiH = np.log10(PhiH)
        hden_Rmax = calculate_hden_Rmax(args.Rmax, n0) # first guess for outer density
        extras.x0 = args.initialx0                     # Fix x0 while we converge Rsonic
        converged = False


        extras.it = 0
        scipy.optimize.newton(find_Rsonic_for_hden, 
                              hden_Rmax, tol=args.tolerance**2, args=(args, extras))
        print "Angle %i converged Rsonic after %i iteration%s" \
            % (theta, extras.it, "s" if extras.it > 1 else "")

        scipy.optimize.newton(find_Rif_for_x0, 
                              args.initialx0, tol=args.tolerance**2, args=(args, extras))
        print "Angle %i converged Rif after %i iteration%s" \
            % (theta, extras.it, "s" if extras.it > 1 else "")


    if args.multiprocess:
        # Clean up everything
        pool.close()
        pool.join()


