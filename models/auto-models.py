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
    R = cmdargs.Rmax - m.ovr.depth/cmdargs.r0
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
    R = cmdargs.Rmax - m.ovr.depth/cmdargs.r0
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

def setup_and_run_cloudy_model(cmdargs, extras, monitor):
    cloudy_cmd = ["time", cmdargs.cloudyexec] 
    # files are simply named after iteration number and density
    filename = "it%.2in%.3ex0%.2f" % (monitor.it, extras.hden, extras.x0)
    print "Calculating " + os.path.join(extras.datapath, filename)
    m = Model(filename, dryrun=False, indir=extras.datapath, outdir=extras.datapath, 
              verbose=True, cloudy_cmd=cloudy_cmd)
    if cmdargs.fastcloudy:
        dust = "FastOrion"
    else:
        dust = "Orion"
    m.write(physical.proplyd(cmdargs.r0, extras.hden, Rmax=cmdargs.Rmax, 
                             W=cmdargs.W, x0=extras.x0, dust=dust))
    m.write(misc.optimize())
    m.write(misc.stopping())
    if cmdargs.fastcloudy:
        m.write(misc.iterate(iterations=0))
    else:
        m.write(misc.iterate())
    m.write(incident.background())
    m.write(incident.star(extras.logPhiH, cmdargs.Tstar, atmosphere=cmdargs.atmosphere, log_g=4.1)) # th1C
    m.write(save.default())
    if cmdargs.savecont:
        m.write("""\
save every last ionizing continuum 1.0 0 ".icont"
""")
    m.run()
    return m
    

class BracketedException(Exception):
    pass

class Monitor(object):
    """
    Class for storing variables that monitor the state of the root-finding process

    For instance, iteration number, list of previous x and y values, etc
    
    None of the variables stored here should materially affect the model calculations
    """
    def __init__(self, startiter=0, onlybracket=False):
        self.it = startiter
        # Dictionary for storing function evaluations with different hden or x0
        self.cache = dict()
        # x-values that most closely bracket the root
        self.xpos = None
        self.xneg = None
        # In the first stage we only want to bracket the root
        self.onlybracket = onlybracket
        

    def update(self, x, y):
        # Update the cache
        self.cache[x] = y
        # Update the brackets
        if y > 0.0:
            if self.xpos is None or y < self.cache[self.xpos]:
                self.xpos = x
        else:
            if self.xneg is None or y > self.cache[self.xneg]:
                self.xneg = x

        print "Current brackets: ", self.xneg, self.xpos
        # Return True if we have brackets
        return not (self.xpos is None or self.xneg is None)
    




def find_Rsonic_for_hden(hden_Rmax, cmdargs, extras, monitor):
    """
    Runs a Cloudy model to find the dimensionless radius R of the sonic point
    
    Principal argument: hden_Rmax

    Other vars are taken from cmdargs (command line parameters) and extras

    it - problematic when called as part of root finding routine
    """
    extras.hden = hden_Rmax
    if hden_Rmax in monitor.cache:
        result = monitor.cache[hden_Rmax]
        print "Using cached result for hden = %.3e" % (hden_Rmax)
    else:
        print "Running cloudy with hden = %.3e" % (hden_Rmax)
        monitor.it += 1
        m = setup_and_run_cloudy_model(cmdargs, extras, monitor)
        result = find_max_sound(m) - 1.0
        bracketed = monitor.update(hden_Rmax, result)
        if monitor.onlybracket and bracketed:
            # Option to give up on looking for this root by the first method
            raise BracketedException
    return result

def find_Rif_for_x0(x0, cmdargs, extras, monitor):
    """
    Runs a Cloudy model to find the dimensionless radius R of the i-front
    
    Principal argument: x0

    Other vars are taken from cmdargs (command line parameters) and extras

    it - problematic when called as part of root finding routine
    """
    extras.x0 = x0
    dR = (0.5*cmdargs.W*extras.x0) / \
        (cmdargs.r0*SIGMA0*extras.hden*cmdargs.Rmax*cmdargs.Rmax*velocity_law_U_of_R(cmdargs.Rmax))
    if x0 in monitor.cache:
        result = monitor.cache[x0]
        print "Using cached result for x0 = %.2f, dR = %.5f" % (x0, dR)
    else:
        print "Running cloudy with x0 = %.2f, dR = %.5f" % (x0, dR)
        monitor.it += 1
        m = setup_and_run_cloudy_model(cmdargs, extras, monitor)
        Rif = find_ifront(m)
        result = Rif - (1.0 - dR)
        bracketed = monitor.update(x0, result)
        if monitor.onlybracket and bracketed:
            # Option to give up on looking for this root by the first method
            raise BracketedException
        print "I-front radius = %.4f, target = %.4f" % (Rif, (1.0 - dR))

    return result





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
    parser.add_argument("--initialhden", type=float, default=0.0,
                        help="Initial guess at outer density (if 0.0 then this will be estimated from Stroemgren condition)")
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

    
    cmdargs = parser.parse_args()

    
    
    if cmdargs.multiprocess:
        raise NotImplementedError # Delete this when it is implemented properly
        pool = multiprocessing.Pool() # by default makes one process per core

    # ID for this proplyd model
    modelid = "%s%.6i-phi%.2f-r%.2f" % (cmdargs.atmosphere, cmdargs.Tstar, 
                                        cmdargs.logPhiH, np.log10(cmdargs.r0))
    # the ID may have to be refined by adding more parameters later

    thetas = np.linspace(0.0, 90.0, cmdargs.ntheta)
    for theta in thetas:
        extras = ExtraPars() # container for miscellaneous parameters
        thetaid = "th%.2i" % (theta)
        # All the files for a given angle now go in the same directory...
        extras.datapath = os.path.join(cmdargs.rootdir, modelid, thetaid)
        PhiH = (10**cmdargs.logPhiH) * (np.cos(theta * DEGREES) + cmdargs.diffuseBeta)
        if theta == 0.0:
            PhiH0 = PhiH
            if cmdargs.initialhden == 0.0:
                # Only estimate the outer density if it was not specified on the command line
                n0 = calculate_nominal_density(PhiH, cmdargs.r0)
                hden_Rmax = calculate_hden_Rmax(cmdargs.Rmax, n0) # first guess for outer density
            else:
                hden_Rmax = cmdargs.initialhden
            hden_Rmax0 = hden_Rmax
        else:
            hden_Rmax = hden_Rmax0 * np.sqrt(PhiH/PhiH0)
            
        extras.logPhiH = np.log10(PhiH)
        extras.x0 = cmdargs.initialx0                     # Fix x0 while we converge Rsonic

        hdenmonitor = Monitor(startiter=0, onlybracket=True)
        try: 
            scipy.optimize.newton(find_Rsonic_for_hden, 
                                  hden_Rmax, tol=cmdargs.tolerance**2, 
                                  args=(cmdargs, extras, hdenmonitor))
        except BracketedException as error:
            print error
            hdenmonitor.onlybracket = False
            scipy.optimize.brentq(find_Rsonic_for_hden, 
                                  a=hdenmonitor.xneg, b=hdenmonitor.xpos, 
                                  xtol=cmdargs.tolerance*hdenmonitor.xneg,
                                  args=(cmdargs, extras, hdenmonitor))
        except RuntimeError as error: 
            print "Ignoring RuntimeError: ", error

        print "Angle %i converged Rsonic after %i iteration%s" \
            % (theta, hdenmonitor.it, "s" if hdenmonitor.it > 1 else "")

        x0monitor = Monitor(startiter=hdenmonitor.it, onlybracket=True)
        try:
            scipy.optimize.newton(find_Rif_for_x0, 
                                  cmdargs.initialx0, tol=cmdargs.tolerance**2, 
                                  args=(cmdargs, extras, x0monitor))
        except BracketedException as error:
            print error
            x0monitor.onlybracket = False
            scipy.optimize.brentq(find_Rif_for_x0,
                                  a=x0monitor.xneg, b=x0monitor.xpos, 
                                  xtol=cmdargs.tolerance*x0monitor.xneg,
                                  args=(cmdargs, extras, x0monitor))
        except RuntimeError as error: 
            print "Ignoring RuntimeError: ", error

        print "Angle %i converged Rif after %i iteration%s" \
            % (theta, x0monitor.it, "s" if x0monitor.it > 1 else "")


    if cmdargs.multiprocess:
        # Clean up everything
        pool.close()
        pool.join()


