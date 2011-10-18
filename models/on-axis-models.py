"""
Run a series of cloudy proplyd models to test the density law


Original version based on ~/Work/Bowshocks/LLobjects/Cloudy/ll1-thon-sd.py

"""

import numpy as np
import scipy
from scipy.optimize import fsolve, bisect
import subprocess
import sys
import multiprocessing

from cloudycontrol.model import Model
from cloudycontrol import incident, save, misc, physical


scriptname = sys.argv[0]

MULTIPROCESS=True           # whether or not to use multiprocess.Pool


rad_list = [
    "WM",
    # "TL",
    ]

# Proplyd parameters for LV2
r0 = 8e14
Rmax = 9.0

def calculate_phi(r0, n0, omega=0.121049705173, alpha=2.6e-13):
    """
    Ionizing flux required for a nominal n0 and r0: Phi = omega r0 n0**2 alpha
    """
    return omega*alpha*r0*n0**2

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
    

# These are the values that Nahiely has been using
hden_Rmax = 1.e4
log_PhiH = 14.54
x0 = 10.8
A = 174.0


# By hand adjustment to try and fix i-front
beta_index = -4.5
hden_Rmax *= 1.0121**beta_index

Tstar_range = [
    3.8e4, 
    # 3.9e4, 
    # 4.0e4,
    ]


cloudy_cmd = ["time", "/Users/will/Work/CLOUDY/SVN/trunk/source/cloudy.exe"] 

if __name__ == '__main__':
    if MULTIPROCESS:
        pool = multiprocessing.Pool() # by default makes one process per core

    for Tstar in Tstar_range:
        for radiation in rad_list:
            filename = "proplyd-test-%s%.2i-phi%.2f-n%.2e-A%.3i" % (radiation, Tstar/1000,
                                                                    log_PhiH, hden_Rmax,
                                                                    A)
            print "Calculating " + filename
            m = Model(filename, verbose=True, cloudy_cmd=cloudy_cmd)
            m.write(physical.proplyd(r0, hden_Rmax, Rmax=Rmax, A=A, x0=x0))
            m.write(misc.optimize())
            m.write(misc.stopping())
            m.write(misc.iterate())
            m.write(incident.background())
            m.write(incident.star(log_PhiH, Tstar, atmosphere=radiation, log_g=4.1)) # th1C
            m.write(save.default())
            # Do full calculation
            if MULTIPROCESS:
                pool.apply_async(m)
            else:
                m.run()

    if MULTIPROCESS:
        pool.close()
        pool.join()


