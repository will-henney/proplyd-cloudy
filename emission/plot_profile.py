import numpy as np
import pyfits
import matplotlib.pyplot as plt
import argparse
import sys
import itertools 

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="""
Plot line profiles from 1D FITS files that were created by extract_aperture.py
    """)

parser.add_argument(
    "lines", type=str, nargs="+",  
    help='Cloudy IDs of emission lines to plot (with _ instead of space)')

parser.add_argument(
    "--raw", action="store_true",
    help='Use the raw fluxes (default is to normalize to peak)')

parser.add_argument(
    "--logscale", action="store_true",
    help="Use logarithmic scale for the fluxes (default is linear)")

parser.add_argument(
    "--vrange", nargs=2, type=float, metavar=("vmin", "vmax"),
    help="Velocity range for plot (default is limits of data)")

parser.add_argument(
    "--legendpos", type=str, default="best",
    help="Manually adjust placement of graph legend")

parser.add_argument(
    "--output", type=str, default="plot_profile.png",
    help="Name of output file (type determined frm suffix)")

parser.add_argument(
    "--fluxmin", type=float, 
    help="Minimum flux for graph (useful for log plots)")

parser.add_argument(
    "--jitter", type=float, default=0.5, 
    help="Amplitude of random jitter to add to velocities (to visually separate the histograms)")

cmd_args = parser.parse_args()

fig = plt.figure()
ax = fig.add_subplot(111)

# This trick is from Avaris on stackoverflow:
# http://stackoverflow.com/questions/7799156/can-i-cycle-through-line-styles-in-matplotlib
linecycler = itertools.cycle(["-","--"])
markcycler = itertools.cycle(["o", "v", "s"])
lwcycler = itertools.cycle([4, 2])
alphacycler = itertools.cycle([0.4, 0.8])

ax.set_color_cycle(
    ["#ff0000", "#00aa00", "#0000ff", "#00ccdd", "#ff00ff", "#ffaa00", "#000000"]
    )

for line in cmd_args.lines:
    fitsfile = "profile-{}.fits".format(line)
    hdu, = pyfits.open(fitsfile)
    # set up U array from WCS header keywords
    k0, u0, du, nu = [hdu.header[kwd] for kwd in "CRPIX1", "CRVAL1", "CDELT1", "NAXIS1"]
    U = u0 + du*(np.arange(nu) - k0) # velocity array
    U += cmd_args.jitter * np.random.uniform() * du
    F = hdu.data                     # flux array
    if not cmd_args.raw:
        F /= F.max()            # optionally scale each flux by its maximum
    # choose log or linear flux axis
    plot = ax.semilogy if cmd_args.logscale else ax.plot
    # add to the plots
    p, = plot(U, F, linestyle='-', label=line, 
              # marker=next(markcycler), 
              lw=next(lwcycler), alpha=next(alphacycler), 
              # drawstyle="steps-mid"
              )
    ax.fill_between(U, F, alpha=0.4, color=p.get_color())

figfile = cmd_args.output
ax.set_xlabel("Velocity, km/s")
ax.set_ylabel("Raw flux" if cmd_args.raw else "Normalized flux")
ax.set_xlim(cmd_args.vrange)
if cmd_args.fluxmin:
    ax.set_ylim(bottom=cmd_args.fluxmin)
# Make a legend with showing the line styles
leg = ax.legend(loc=cmd_args.legendpos, prop=dict(size=10))
fig.savefig(figfile)
