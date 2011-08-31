"""
Plot various physical variables across the ionization front

Use a split double-log scale around the point where x = 0.5

Based on an earlier incarnation in proplyd.org
"""

import numpy, sys, os
sys.path.append("../src"); import claudia
import matplotlib.pyplot as plt
import argparse


# Specify the Cloudy model on the command line
parser = argparse.ArgumentParser(description="Plot physical variables across the ionization front")
parser.add_argument("--model", "-m", type=str, 
                    help="Name of Cloudy model")
args = parser.parse_args()              # we can now use args.runid, args.itime


  
r0 = 8.e14
Rmax = 9.0
deltascale = 1.e11

# Parameters for ionization transition, specified in .in files
x0 = 2.3                                  # R = 1 occurs here
A = 15.0

# photoionization cross section at 1 Rydberg
sigma = 6.3e-18

  
plt.figure(1)

m = claudia.CloudyModel(args.model, indir="in", outdir="out")
z = m.ovr.depth
hden = 10**m.ovr.hden
eden = 10**m.ovr.eden
efrac = eden/hden
ifrac = 10**m.ovr.HII
T = 10**m.ovr.Te
R = Rmax - z/r0
sound = m.pre.cadwind_kms

dz = z[1:] - z[:-1]
dz = numpy.array(dz.tolist() + [dz[-1]])
colden = numpy.cumsum(hden*dz)

Tmax = T.max()
i1 = T.argmax()
i2 = len(R[R>=1.0])
n0 = hden[i2]
U = n0 / (R**2 * hden)

i5 = len(R[ifrac>=0.5])                 # position of ifrac = 0.5

# calculate the dimensionless ionization thickness in the approximate model
rho_Rmax = hden[0]
dR = A/(r0*sigma*rho_Rmax*Rmax*Rmax*U[0])

print "-"*72
print "Tmax at R = %.4f, ifrac = %.4f" % (R[i1], ifrac[i1])
print "@ R = 1.0: T/Tmax = %.4f, ifrac = %.4f, U = %.4f" % (T[i2]/Tmax, ifrac[i2], U[i2])
print "50%% ionization at R = %.4f" % (R[i5])

delta = (z[i5] - z)/deltascale
plt.plot(delta, U, 'k-', label='velocity')                  # cloudy version
plt.plot(delta, T/1.e4, 'r-', label='T / 10^4 K')
plt.plot(delta, eden/n0, 'k--', label='n_e / n_0')
plt.plot(delta, sound/sound[i2], 'r--', label='c / c_m')
plt.plot(delta, efrac, 'b-', label='n_e / n_H')
plt.plot(delta, 5e-22*colden, 'b--', label='A_V')

plt.xlabel('[r - r(x=0.5)] / %.0e cm' % (deltascale))
plt.ylabel('')
plt.title('Proplyd ionization front structure')
plt.grid(True)
plt.xscale('symlog')
plt.axis([-0.3*r0/deltascale, Rmax*r0/deltascale, 0.0, 2.8])
plt.legend(loc="upper left")
plt.savefig("ifront-structure-%s.png" % (args.model))
