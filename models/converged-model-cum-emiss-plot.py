import numpy, sys, os
sys.path.append("/Users/will/Work/Nahiely/proplyd-cloudy/src")
import claudia
import matplotlib.pyplot as plt
  


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
    sound = numpy.sqrt(3./5.)*m.pre.cadwind_kms

    # Now get the delta z's directly from the ".dr" file
    # We need to skip the values marked as zone 0 since they are spurious
    dz = m.dr.dr[m.dr.zone != 0] 
    assert len(dz) == len(z)    # check we got it right
    colden = numpy.cumsum(hden*dz)

    linelist = [
        'O__3__5007A',
        'TOTL__4363A',
        'C__3__1910A',
        'H__1__6563A',
        'N__2__6584A',
        'O_II__3726A',
        'S_II__6731A',
        'O__1__6300A',
        ]

    emlines = dict()
    emsum = dict()
    for id in linelist:
        emissivity = 10**m.em[id]
        cumemiss = numpy.cumsum(R**2 * emissivity*dz)
        emsum[id] = cumemiss.max()
        emlines[id] = cumemiss/cumemiss.max()

    Tmax = T.max()
    i1 = T.argmax()
    i2 = len(R[R>=1.0])
    n0 = hden[i2]
    U = n0 / (R**2 * hden)

    # Heavy element ion fracs
    oplus = m.ion_o.O__2
    nplus = m.ion_n.N__2
    splus = m.ion_s.S__2

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

    print "-"*72
    
    delta = (z[i5] - z)/deltascale

    v = U*sound
    centiles = [10, 25, 50, 75, 90]
    for emid in linelist:
        radii = [R[emlines[emid] > 0.01*cent].max() for cent in centiles]
        formatted_radii = ["R(%.2i)=%5.3f" % (cent, rad) for cent, rad in zip(centiles, radii)]
        print "%s:" % (emid), "Flux = %.2e," % (emsum[emid]), ", ".join(formatted_radii)
        velocities = [v[emlines[emid] > 0.01*cent].max() for cent in centiles]
        formatted_velocities = ["V(%.2i)=%5.2f" % (cent, vel) for cent, vel in zip(centiles, velocities)]
        print "           :", "  F/Ha = %.4f," % (emsum[emid]/emsum['H__1__6563A']), ", ".join(formatted_velocities)

    
    for emid in linelist:
        plt.plot(delta, emlines[emid], label=emid)

    plt.xlabel('[r - r(x=0.5)] / %.0e cm' % (deltascale))
    plt.ylabel('')
    plt.title('Cumulative emission line flux')
    plt.grid(True)
    plt.xscale('symlog')
    plt.axis([-0.3*r0/deltascale, Rmax*r0/deltascale, 0.0, 1.1])
    plt.legend(loc="lower left")
    plt.savefig("emissplot-%s.png" % (modelid))

if __name__ == '__main__':
    import warnings, argparse
    warnings.filterwarnings("ignore")

    parser = argparse.ArgumentParser(
        description="Make emissivity plots of Cloudy model i-front structure on split bi-log scale",
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
    args = parser.parse_args()
    

    r0, deltascale = args.r0, args.scale
    A, x0 = args.A, args.x0
    plot_vars(args.id)
