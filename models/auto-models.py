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


scriptname = sys.argv[0]

MULTIPROCESS=False           # whether or not to use multiprocess.Pool

DEGREES = np.pi/180.0           # conversion from degrees to radians

# Proplyd parameters for LV2
r0 = 8e14
Rmax = 9.0
beta = 0.004                     # Diffuse/stellar flux ratio

def calculate_phi(r0, n0, omega=0.121049705173, alpha=2.6e-13):
    """
    Ionizing flux required for a nominal n0 and r0: Phi = omega r0 n0**2 alpha
    """
    return omega*alpha*r0*n0**2

def calculate_nominal_density(phi, r0, omega=0.121049705173, alpha=2.6e-13):
    """
    Complementary calculation of density as a function of ionizing flux
    """
    return np.sqrt(phi / (omega*alpha*r0))

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
    

def find_max_sound(model):
    "Find the radius at which sound speed is maximum"
    R = 1.0                     # Obviously, needs improvement!
    return R


# theta is the angle of incidence of the stellar radiation at the i-front
thetas = [0.0, 30.0, 60.0, 90.0]

# These are the values that Nahiely has been using
hden_Rmax = 1.e4
PhiH_star = 10**14.54
x0 = 10.8
A = 174.0



# Automatic adjustment to try and fix i-front
beta_index = -4.5
# hden_Rmax *= 1.0121**beta_index

Tstar = 3.8e4
model_atmosphere = "WM"


cloudy_cmd = ["time", "/Users/will/Work/CLOUDY/SVN/trunk/source/cloudy.exe"] 

rootdir = "proplyd-auto-models"

tolerance = 1.e-3

if __name__ == '__main__':
    if MULTIPROCESS:
        pool = multiprocessing.Pool() # by default makes one process per core

    # ID for this proplyd model
    modelid = "%s%.2i-phi%.2f" % (model_atmosphere, Tstar/1000, np.log10(PhiH_star))
    # the ID may have to be refined by adding more parameters later

    for theta in thetas:
        thetaid = "th%.2i" % (theta)
        PhiH = PhiH_star*(np.cos(theta * DEGREES) + beta)
        n0 = calculate_nominal_density(PhiH, r0)
        log_PhiH = np.log10(PhiH)
        hden_Rmax = calculate_hden_Rmax(Rmax, n0) # first guess for outer density

        converged = False
        it = 0
        while not converged:
            it += 1
            # All the files for a given angle now go in the same directory...
            datapath = os.path.join(rootdir, modelid, thetaid)
            # files are simply named after iteration number and density
            filename = "it%.2in%.3e" % (it, hden_Rmax)
            print "Calculating " + os.path.join(datapath, filename)

            m = Model(filename, dryrun=True, indir=datapath, outdir=datapath, 
                      verbose=True, cloudy_cmd=cloudy_cmd)
            m.write(physical.proplyd(r0, hden_Rmax, Rmax=Rmax, A=A, x0=x0))
            m.write(misc.optimize())
            m.write(misc.stopping())
            m.write(misc.iterate())
            m.write(incident.background())
            m.write(incident.star(log_PhiH, Tstar, atmosphere=model_atmosphere, log_g=4.1)) # th1C
            m.write(save.default())
            # Do full calculation
            if MULTIPROCESS:
                pool.apply_async(m)
            else:
                m.run()

            # CHECK IF CONVERGED
            Rsonic = find_max_sound(m)
            if abs(Rsonic - 1.0) > tolerance:
                print "Sonic point is in wrong place (R=%.2f) - adjust density and try again"
                # IF NOT, FIND CORRECTED OUTER DENSITY
                hden_Rmax *= Rsonic**beta_index
            else:
                converged = True

        print "Angle %i converged in %i iteration%s" % (theta, it, "s" if it > 1 else "")

    if MULTIPROCESS:
        pool.close()
        pool.join()


