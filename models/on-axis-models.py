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

scriptname = sys.argv[0]

MULTIPROCESS=True           # whether or not to use multiprocess.Pool

cloudy_input = \
    '''* Cloudy input file written by %(scriptname)s
title Proplyd on-axis test (atmos=%(radiation)s, Tstar=%(Tstar).1f, phi(H) = %(log_PhiH).2f, den(Rmax)=%(hden_Rmax).1f, A=%(A).1f)
set save prefix "%(filename)s"
dlaw 53, %(r0).1e %(Rmax).1f %(hden_Rmax).1f %(A).1f %(x0).2f
iterate
stop temperature 5 linear
stop efrac = 0.01
'''

optimize_input = \
'''no molecules
// DISABLED: no level 2 lines
'''

cloudy_savefiles_input = \
'''* Output options
set save hash "return"
set save flush
save last overview ".ovr"
save last physical conditions ".phy"
save last grain charge ".grc"
save last grain drift velocity ".grv"
save last grain temperature ".grt"
save last PDR ".pdr"
save last element oxygen ".ion_O"
save last element nitrogen ".ion_N"
save last element sulphur ".ion_S"
save last element carbon ".ion_C"
save last dr ".dr"
save last continuum ".cont"
save last pressure ".pre"
save last heat ".heat" 
save last cool ".cool" 
save last lines, emissivity ".em"
H  1  4861A 
H  1  6563A 
TOTr  5199A 
N  2  6584A 
O  1  6300A 
O II  3726A 
TOTL  4363A 
O  3  4959A 
O  3  5007A 
S II  6731A 
S II  6716A 
S  3  6312A
He 1  5876A
Ne 2 12.81m
Ne 3 15.55m
O  2 4651A
C  3 1910A
C  3 1907A
end of lines
'''

baseline_radiation_input = """* Stuff that is always there
cosmic ray background 
cmb
table ism // Diffuse interstellar field
"""

orion_dust_abundances_input = """* Orion nebula abundances plus Orion dust
abundances H II region no grains
grains Orion
grains PAH
set PAH "H" // Only have PAH in the neutral gas
"""

cloudy_cmd = ["time", "/Users/will/Work/CLOUDY/SVN/trunk/source/cloudy.exe"] 


blackbody_radiation_input = """* Photoionization equilibrium
black body, T=%(Tstar).1f K 
""" 
wmbasic_radiation_input = """* WMBasic (Pauldrach et al. 2001) non-LTE, line-blanketed, and wind-blanketed hot stars
table star wmbasic %(Tstar).1f %(log_g).1f %(log_Z).1f 
"""

atlas_radiation_input = """* Atlas (Castelli & Kurucz 2004) LTE, plane-parallel, hydrostatic model atmospheres
table star atlas odfnew Z+%(log_Z).1f %(Tstar).1f %(log_g).1f
"""

tlusty_radiation_input = """* Tlusty (Lanz & Hubeny 2003) non-LTE, line-blanketed, plane-parallel, hydrostatic O and B star SEDs.
table star tlusty OBstar 3-dim %(Tstar).1f %(log_g).1f %(log_Z).1f
"""

kurucz_radiation_input = """* Kurucz (1979) similar to Atlas but obsolete - used only for comparison with Baldwin (1991)
table star kurucz %(Tstar).1f
"""


spectrum_shape = dict(
    BB = blackbody_radiation_input,
    WM = wmbasic_radiation_input,
    AT = atlas_radiation_input,
    TL = tlusty_radiation_input,
    KU = kurucz_radiation_input
    )

rad_list = [
    "WM",
    # "TL",
    ]

def stellar_spectrum(log_phiH, Tstar, log_g=4.0, log_Z=0.0, model='BB'):
    s = spectrum_shape[model] % dict(Tstar=Tstar, log_g=log_g, log_Z=log_Z)
    s += "phi(H) %.2f\n" % (log_phiH)
    return s


abundances_input = orion_dust_abundances_input

def run_cloudy(filename):
    infile = file('in/%s.in' % filename, 'r')
    outfile = file('out/%s.out' % filename, 'w')
    subprocess.Popen(cloudy_cmd, 
                     stdin=infile, 
                     stdout=outfile,
                     stderr=subprocess.STDOUT,
                     cwd='out/'
                     ).wait()
    if MULTIPROCESS:
        print multiprocessing.current_process().name, " finished with ", filename


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



if __name__ == '__main__':
    if MULTIPROCESS:
        pool = multiprocessing.Pool() # by default makes one process per core

    for Tstar in Tstar_range:
        for radiation in rad_list:
            filename = "proplyd-test-%s%.2i-phi%.2f-n%.2e-A%.3i" % (radiation, Tstar/1000,
                                                                    log_PhiH, hden_Rmax,
                                                                    A)
            print "Calculating " + filename
            with file('in/%s.in' % filename, 'w') as infile:
                infile.write(cloudy_input % locals())
                infile.write(optimize_input)
                infile.write(abundances_input)
                infile.write(baseline_radiation_input)
                infile.write(stellar_spectrum(log_PhiH, Tstar, 4.1, model=radiation)) # th1C
                infile.write(cloudy_savefiles_input)
            # Do full calculation
            if MULTIPROCESS:
                pool.apply_async(run_cloudy, (filename,))
            else:
                run_cloudy(filename)

    if MULTIPROCESS:
        pool.close()
        pool.join()


