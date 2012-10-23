import numpy, sys, os
# ugly hack to make sure the claudia src folder is in sys.path
sys.path.append(
    os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 
        "src"
        )
    )
import claudia
import matplotlib.pyplot as plt
import matplotlib 

matplotlib.rcParams['lines.linewidth'] = 3


# Parameters for ionization transition, specified in .in files
x0 = 10.8                                  # R = 1 occurs here
A = 174.0
Rmax = 9.0

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
    nfrac = 10**m.ovr.HI
    T = 10**m.ovr.Te
    R = Rmax - z/r0
    sound = m.pre.cadwind_kms

    dz = z[1:] - z[:-1]
    dz = numpy.array(dz.tolist() + [dz[-1]])
    colden = numpy.cumsum(hden*dz)

    Tmax = T.max()
    i1 = T.argmax()
    i2 = len(R[R>=1.0])
    print i2, len(hden)
    n0 = hden[i2]
    U = n0 / (R**2 * hden)

    # Heavy element ion fracs
    oplus = m.ion_o.O__2
    nplus = m.ion_n.N__2
    splus = m.ion_s.S__2
    heplus = 10**m.ovr.HeII

    # Calculate the optical depths

    # First at the Lyman limit
    sig0 = 6.e-18 
    tau0 = sig0*numpy.cumsum(hden*nfrac*dz)

    # Second, the "mean" ionizing optical depth
    # TODO

    # Heating/cooling rates
    # We want the values per n_e n_i in units of 1.e-24 erg cm^3 / s
    heat24 = m.phy.Htot/(1.e-24*eden*hden*ifrac) 

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
    plt.plot(delta, sound/sound[i2], 'r-.', label='c / c_m')
    plt.plot(delta, eden/n0, 'k--', label='n_e/n_0')
    plt.plot(delta, heat24/3.0, 'r--', label='Heat / 3.e-24')
    plt.plot(delta, efrac, 'b-', label='n_e/n_H')
    plt.plot(delta, tau0/10.0, 'b--', label='tau0/10')

    plt.plot(delta, oplus, 'g-', label='O+/O')
    plt.plot(delta, nplus, 'g--', label='N+/N')
    plt.plot(delta, heplus, 'y-', label='He+/He')
    plt.plot(delta, splus, 'y--', label='S+/S')

    for l in plt.gca().lines:
        l.set_alpha(.7)

    plt.xlabel('[r - r(x=0.5)] / %.0e cm' % (deltascale))
    plt.ylabel('')
    plt.title('Proplyd ionization front structure')
    plt.grid(True)
    plt.xscale('symlog')
    plt.axis([-0.3*r0/deltascale, Rmax*r0/deltascale, 0.0, 2.8])
    plt.legend(loc="upper left", prop={'size':10})
    plt.savefig("thermoplot-%s.png" % (modelid), dpi=args.dpi)

if __name__ == '__main__':
    import warnings, argparse
    warnings.filterwarnings("ignore")

    parser = argparse.ArgumentParser(
        description="Make thermal plots of Cloudy model i-front structure on split bi-log scale",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument("id", type=str,
                        help='Name of Cloudy model')
    parser.add_argument("--x0", type=float, default=10.8,
                        help='I-front x0 parameter')
    parser.add_argument("--A", type=float, default=174.0,
                        help='I-front A parameter')
    parser.add_argument("--r0", type=float, default=8e14,
                        help='I-front radius')
    parser.add_argument("--scale", type=float, default=1.e11,
                        help='Length scale for plot (should be roughly r0/1.e4)')
    parser.add_argument("--dpi", type=int, default=300,
                        help='Resolution of PNG file (values >= 600 give multi-MB files)')
    args = parser.parse_args()
    

    r0, deltascale = args.r0, args.scale
    A, x0 = args.A, args.x0
    plot_vars(args.id)
