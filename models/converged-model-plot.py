import numpy, sys, os
sys.path.append("/Users/will/Work/Nahiely/proplyd-cloudy/src")
import claudia
import matplotlib.pyplot as plt
  
r0 = 8.e14
Rmax = 9.0
deltascale = 1.e11


# Parameters for ionization transition, specified in .in files
x0 = 10.8                                  # R = 1 occurs here
A = 174.0

# photoionization cross section at 1 Rydberg
sigma = 6.3e-18

  
def ionfrac(R, dR):
    "Smooth change of ion fraction at the i-front"
    x = x0*(R - 1.0 + dR)/dR
    return 0.5*(numpy.tanh(x) + 1)

def ionfrac2(R, f1=0.0, f2=1.0):
    "Smooth change of ion fraction between f1 on neutral side and f2 on ionized side"
    x = x0*(R - 1.0 + dR)/dR
    return f1 + (f2 - f1)*0.5*(numpy.tanh(x) + 1)

def ioniz(x):
    "Degreee pf ionization: log x / (1 - x)"
    return numpy.log10(x/(1.-x))

def plot_vars(modelid):

    plt.figure(1)
    m = claudia.CloudyModel(modelid)
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

    # Approximate ionization fraction that we use in dense_fabden.cpp
    xapprox = ionfrac(R, dR)

    print "-"*72
    print "model %s: Tmax at R = %.4f, ifrac = %.4f" % (modelid, R[i1], ifrac[i1])
    print "@ R = 1.0: T/Tmax = %.4f, ifrac = %.4f, U = %.4f" % (T[i2]/Tmax, ifrac[i2], U[i2])
    print "I-front thickness in approx model: dR = %.4f" % (dR)

    delta = (z[i5] - z)/deltascale
    plt.plot(delta, U, 'k-', label='velocity')                  # cloudy version
    plt.plot(delta, T/1.e4, 'r-', label='T / 10^4 K')
    plt.plot(delta, eden/n0, 'k--', label='n_e/n_0')
    plt.plot(delta, sound/sound[i2], 'r--', label='c / c_m')
    plt.plot(delta, efrac, 'b-', label='n_e/n_H')
    plt.plot(delta, 5e-22*colden, 'b--', label='A_V')
    plt.plot(delta, xapprox, '.', label='approx n+/n_H')

    plt.xlabel('[r - r(x=0.5)] / %.0e cm' % (deltascale))
    plt.ylabel('')
    plt.title('Proplyd ionization front structure')
    plt.grid(True)
    plt.xscale('symlog')
    plt.axis([-0.3*r0/deltascale, Rmax*r0/deltascale, 0.0, 2.8])
    plt.legend(loc="upper left")
    plt.savefig("varplot-%s.png" % (modelid))

if __name__ == '__main__':
    import warnings, argparse
    warnings.filterwarnings("ignore")

    parser = argparse.ArgumentParser(
        description="Make diagnostic plots of Cloudy model i-front structure on twin-log scale",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument("id", type=str,
                        help='Name of Cloudy model')
    parser.add_argument("--x0", type=float, default=10.8,
                        help='I-front x0 parameter')
    parser.add_argument("--A", type=float, default=174.0,
                        help='I-front A parameter')
    args = parser.parse_args()
    
    A, x0 = args.A, args.x0
    plot_vars(args.id)
