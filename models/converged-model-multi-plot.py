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


params = {
    "lines.linewidth": 3,
    "font.family": "serif",
    "text.usetex": True,
    "text.latex.preamble": [r"\usepackage[varg]{txfonts}"],
    "figure.figsize": (5, 10),
    }
matplotlib.rcParams.update(params)


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
    # convert adiabatic -> isothermal sound speed
    sound = m.pre.cadwind_kms / numpy.sqrt(5./3.)

    dz = z[1:] - z[:-1]
    dz = numpy.array(dz.tolist() + [dz[-1]])
    colden = numpy.cumsum(hden*dz)

    Tmax = T.max()
    i1 = T.argmax()
    i2 = len(R[R>=1.0])
    i3 = sound.argmax()
    print i2, len(hden)
    n0 = hden[i2]
    U = n0 / (R**2 * hden)

    # Heavy element ion fracs
    oplus = m.ion_o.O__2
    nplus = m.ion_n.N__2
    splus = m.ion_s.S__2
    cplus = m.ion_c.C__2
    heplus = 10**m.ovr.HeII
    
    # Try to read in more ions if we can
    try:
        neplus = m.ion_ne.Ne_2
        clplus = m.ion_cl.Cl_2
    except IndexError:
        neplus = None
        clplus = None

    # Calculate the optical depths

    # First at the Lyman limit
    sig0 = 6.e-18 
    tau0 = sig0*numpy.cumsum(hden*nfrac*dz)

    # Visual extinction
    sigV = 5e-22
    tauV = sigV*colden

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
    print "sound max at {:.4f}".format(R[i3])
    print "@ R = 1.0: T/Tmax = %.4f, ifrac = %.4f, U = %.4f" % (T[i2]/Tmax, ifrac[i2], U[i2])
    print "I-front thickness in approx model: dR = %.4f" % (dR)

    if args.symlog:
        delta = (z[i5] - z)/deltascale
    else:
        delta = (z[-1] + deltascale) - z
        print "r - r_min at sonic point is {:.2e}".format(delta[i2])
        print "r - r_min at x=0.5 is {:.2e}".format(delta[i5])
        print "r - r_min at c_max is {:.2e}".format(delta[i3])

    if args.density_scale: 
        nscale = args.density_scale
        exponent = int(numpy.log10(nscale))
        abcissa = nscale/10**exponent
        if int(abcissa) == 1:
            snscale = r"10^{{{:d}}}~\mathrm{{cm^{{-3}}}}".format(exponent)
        else:
            snscale = r"{:.0f} \times 10^{{{:d}}}~\mathrm{{cm^{{-3}}}}".format(abcissa, exponent)
    else:
        nscale = n0
        snscale = "n_0"

    hasEmissivityPanel = not args.emlinefile is None

    numCols = 1
    numRows = 3 if hasEmissivityPanel else 2 

    plt.subplot(numRows, numCols, 1)
    plt.plot(delta, U * sound[i2] / 10., 'k-', label='\(\mathrm{velocity} / 10~\mathrm{km~s^{-1}}\)')                  # cloudy version
    plt.plot(delta, T/1.e4, 'r-', label=r'\(T / 10^4\) K')
    plt.plot(delta, sound/10., 'r--', label=r'\(c / 10~\mathrm{km~s^{-1}}\)')
    plt.plot(delta, eden/nscale, 'k--', label=r'\(n_\mathrm{{e}}/{}\)'.format(snscale))
    plt.plot(delta, hden/nscale, 'm-', label=r'\(n_\mathrm{{H}}/{}\)'.format(snscale))
    if args.heat:
        plt.plot(delta, heat24/3.0, 'r-.', label='Heat / 3.e-24')
    plt.plot(delta, efrac, 'b-', label=r'\(n_\mathrm{e}/n_\mathrm{H}\)')
    plt.plot(delta, tau0/10.0, 'b--', label=r'\(\tau_0/10\)')
    plt.plot(delta, 10*tauV, 'g--', label=r'\(10 \times \tau_V\)')

    plt.ylim(0.0, args.ymax)
    if args.symlog:
        plt.xscale('symlog')
        plt.xlim(-0.3*r0/deltascale, Rmax*r0/deltascale)
    else:
        plt.xscale('log')
        plt.xlim(deltascale, Rmax*r0)
    plt.grid(True)
    for l in plt.gca().lines:
        l.set_alpha(.7)
    plt.legend(loc="upper right", title="Physical variables", ncol=3, prop={'size':9})

    plt.subplot(numRows, numCols, numRows)
    plt.plot(delta, oplus, 'g-', label=r'\(\mathrm{O^+\!/\,O}\)')
    plt.plot(delta, nplus, 'g--', label=r'\(\mathrm{N^+\!/\,N}\)')
    plt.plot(delta, heplus, 'y-', label=r'\(\mathrm{He^+\!/\,He}\)')
    plt.plot(delta, splus, 'y--', label=r'\(\mathrm{S^+\!/\,S}\)')
    plt.plot(delta, cplus, 'm-', label=r'\(\mathrm{C^+\!/\,C}\)')
    if not neplus is None:
        plt.plot(delta, neplus, 'b-', label=r'\(\mathrm{Ne^+\!/\,Ne}\)')
    # if not clplus is None:
    #     plt.plot(delta, clplus, 'r-', label=r'\(\mathrm{Cl^+\!/\,Cl}\)')

    for l in plt.gca().lines:
        l.set_alpha(.7)

    plt.ylabel('')
    plt.ylim(0.0, 1.25)
    if args.print_title:
        plt.title('Proplyd ionization front structure')
    plt.grid(True)
    if args.symlog:
        plt.xlabel('[r - r(x=0.5)] / %.0e cm' % (deltascale))
        plt.xscale('symlog')
        plt.xlim(-0.3*r0/deltascale, Rmax*r0/deltascale)
    else:
        plt.xlabel(r'\((r - r_\mathrm{min})\) / cm')
        plt.xscale('log')
        plt.xlim(deltascale, Rmax*r0)
        

    plt.legend(loc="upper right", ncol=2, title="Ion fractions", prop={'size':9})

    colors = dict(O="g", H="r", N="b")

    if hasEmissivityPanel:
        plt.subplot(numRows, numCols, 2)
        linelist = [lineid for lineid in args.emlinefile.read().split("\n") if lineid]
        nlines = len(linelist)
        for iline, lineid in enumerate(linelist):
            color = colors.get(lineid[0], "k")
            emissivity = 10**m.em[lineid]
            cumemiss = numpy.cumsum(R**2 * emissivity*dz)
            cumemiss /= cumemiss.max()
            centiles = {fraction: delta[len(cumemiss[cumemiss <= fraction])] 
                        for fraction in [0.1, 0.25, 0.5, 0.75, 0.9]}
            plt.plot([centiles[0.1], centiles[0.9]], [iline, iline], '-', lw=2, c=color)
            plt.plot([centiles[0.25], centiles[0.75]], [iline, iline], '-|', lw=4, c=color)
            plt.plot(centiles[0.5], iline, 'o', c=color)
            plt.text(deltascale, iline-0.2, 
                     r"\verb|{}|".format(lineid), 
                     fontdict=dict(size="x-small"))

        plt.xlabel('')
        if args.symlog:
            plt.xscale('symlog')
            plt.xlim(-0.3*r0/deltascale, Rmax*r0/deltascale)
        else:
            plt.xscale('log')
            plt.xlim(deltascale, Rmax*r0)
        plt.ylabel('')
        plt.ylim(-0.5, nlines-0.5)
        plt.axis("off")
        plt.grid(True)

    plt.tight_layout()          # tighten things up
    if args.png:
        plt.savefig("multiplot-%s.png" % (modelid), dpi=args.dpi)
    else:
        plt.savefig("multiplot-%s.pdf" % (modelid))

if __name__ == '__main__':
    import warnings, argparse
    warnings.filterwarnings("ignore")

    parser = argparse.ArgumentParser(
        description="Make multi-panel plots of Cloudy model i-front structure",
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
    parser.add_argument("--density-scale", type=float, default=None,
                        help='Density scale for plot (if None, use n0)')
    parser.add_argument("--ymax", type=float, default=3.5,
                        help='Upper limit of dimensionless y axis')
    parser.add_argument("--dpi", type=int, default=300,
                        help='Resolution of PNG file (values >= 600 give multi-MB files)')
    parser.add_argument("--emlinefile", type=file, default=None,
                        help='''Name of a file that lists emission
                        lines to plot.   If this option is specified,
                        then a candlestick plot of the emissivities
                        will be added to the figure''') 
    parser.add_argument("--symlog", action="store_true",
                        help='Use symmetric log scale on the x-axis')
    parser.add_argument("--heat", action="store_true",
                        help='Plot the heating rate too')
    parser.add_argument("--print-title", action="store_true",
                        help='Print a title above the plot')
    parser.add_argument("--png", action="store_true",
                        help='Write a PNG file instead of a PDF')
    args = parser.parse_args()
    

    r0, deltascale = args.r0, args.scale
    A, x0 = args.A, args.x0
    plot_vars(args.id)
