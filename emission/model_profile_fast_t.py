import numpy as np
import scipy
import sys, os, glob
import argparse
import pyfits
import fastcube
import fastcube_gauss

sys.path.append("../../../src")       # make sure we can find claudia.py
import claudia

# Avoid verbose error messages from numpy during the reading of the Cloudy files
import warnings, argparse
warnings.filterwarnings("ignore")


def parse_command_args():

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""
    Calculate 3-d position-position-velocity cube of rebinned Cloudy proplyd model. 
    """)

    parser.add_argument(
        "--modeldir", type=str, default=".",  
        help='Path to proplyd directory model')
    parser.add_argument(
        "--nphi", type=int, default=400,
        help='Number of angles phi')
    parser.add_argument(
        "--nvel", type=int, default=50,
        help='Number of velocity bins')
    parser.add_argument(
        "--inc", type=float, default=0.0,
        help='Inclination in degrees')

    parser.add_argument(
        "--r0", type=float, default=8.0e14,
        help='Radius of proplyd ionization front in cm')
    parser.add_argument(
        "--Rmax", type=float, default=9.0,
        help='Radius of outer face of cloudy model in units of r0')

    parser.add_argument(
        "--window", type=float, nargs=4, default=(-1.0, -2.0, 3.0, 2.0),
        metavar=("X1", "Y1", "X2", "Y2"),
        help='Ounding box of view window (in units of r0) for the position-position-velocity cube')
    parser.add_argument(
        "--pixelsize", type=float, default=0.1, 
        help='Size of pixels (in units of r0) for the PPV cube')

    parser.add_argument(
        "--rebin", type=str, default="linear-101",
        help='Load rebinned model from "rebin-REBIN"')

    parser.add_argument(
        "--onlylinesin", type=str, default=None,
        help='Only create cubes for lines listed in this file (one per line)')

    parser.add_argument(
        "--gauss", action="store_true",  
        help='Calculate thermal broadening of lines (MAY BE SLOW!)')


    return parser.parse_args()



def find_awt(emline):
    "Find the atomic weight from the Cloudy emission line label"
    symbol_dict = dict(
        H=1.0, HE=4.0, C=12.0, N=14.0, O=16.0, NE=20.2, MG=24.3,
        SI=28.1, S=32.1, CL=35.5, AR=40.0, CA=40.1, FE=55.8, CO=58.9, NI=58.7
        )
    special_dict = dict(TOTL__4363A="O", TOTL__2326A="O", TOTL__5199A="N")
    if emline in special_dict.keys():
        element = special_dict[emline]
    else:
        element = emline.split("_")[0]

    return symbol_dict[element]

        


def read_rebinned_models():
    """
    Read the models that are written with models/rebin-models.py
    """
    def read_fits(varid):
        """
        Read variable VARID from fits file of interpolated model structure
        """
        fitsname = os.path.join(
            cmd_args.modeldir,
            "rebin-%s" % (cmd_args.rebin),
            "%s.fits" % (varid)
            )
        return pyfits.open(fitsname)[0].data

    # Find list of rebinned models that can be used
    rebindirs = [
        os.path.basename(p) for p in 
        glob.glob(os.path.join(cmd_args.modeldir, "rebin-*"))
        ]
    Mu = read_fits("mu")[:,0]
    R = read_fits("R")[0,:]
    i2 = len(R[R<1.0])             # find index where R = 1
    Density2D = read_fits("ovr-hden")
    validmask = Density2D != -1.0
    Density2D = 10**Density2D
    sound_speed2D = np.sqrt(3./5.)*read_fits("pre-cadwind_kms")
    n0 = Density2D[:, i2:i2+1]
    c0 = sound_speed2D[:, i2:i2+1] 
    V = c0 * n0 / (R**2 * Density2D) # V should be 2-dimensional

    T = 10**read_fits("ovr-Te")

    emfiles_found = [
        os.path.basename(p) for p in 
        glob.glob(
            os.path.join(
                cmd_args.modeldir,
                "rebin-%s" % (cmd_args.rebin),
                "em-*.fits"
                )
            )
        ]
    emlines_found = [ s.split("-")[1].split(".")[0] for s in emfiles_found ]

    if cmd_args.onlylinesin:
        # Read list of selected emission lines if asked for
        with open(cmd_args.onlylinesin) as f:
            selectedlines = f.read().split('\n')
        # But make sure that H beta is in the list
        if not "H__1__4861A" in selectedlines:
            selectedlines.append("H__1__4861A")
        print "Selected lines: ", selectedlines
        # Only calculate lines that are not in the selected list
        emlines = list()
        emfiles = list()
        for emline, emfile in zip(emlines_found, emfiles_found):
            if emline in selectedlines:
                emlines.append(emline)
                emfiles.append(emfile)
        print "Lines to be used: ", emlines
    else:
        # Otherwise, use all the lines whose emissivity files can be found
        emlines = emlines_found
        emfiles = emfiles_found
        

    # list of emissivity arrays
    Emissivities = [10**read_fits("em-%s" % (emline)) for emline in emlines]
    # convert in 3D array:
    Ems = np.concatenate( [arr[None,:,:] for arr in Emissivities] )

    # find atomic weights
    awt = [find_awt(emline) for emline in emlines]

    return emlines, Mu, R, V, T, Ems, awt, validmask






def write_datacubes():
    """Write each datacube to a FITS file """
    def add_WCS_keyvals(iaxis, keyvaldict):
        """Add in a whole load of keywords for the same axis"""
        for k, v in keyvaldict.items():
            hdu.header.update(k+str(iaxis), v)
    # write the data cubes as FITS files
    for emline, cube in DataCubes.items():
        fitsfilename = "cube-{}-{}.fits".format(emline, idstring)
        # write datacube to the HDU
        hdu = pyfits.PrimaryHDU(cube)
        # write WCS headers to the HDU
        add_WCS_keyvals(3, dict(crpix=1, crval=umin+0.5*du, cdelt=du, 
                                ctype="VELOCITY", cunit="km/s    ") )
        add_WCS_keyvals(2, dict(crpix=1, crval=ymin+0.5*dy, cdelt=dy, 
                                ctype="YOFFSET ", cunit="r0      ") )
        add_WCS_keyvals(1, dict(crpix=1, crval=xmin+0.5*dx, cdelt=dx, 
                                ctype="XOFFSET ", cunit="r0      ") )
        hdu.writeto(fitsfilename, clobber=True)


def write_lineprofile_table():
    import matplotlib.mlab as mlab
    datatable = np.rec.fromarrays([PerfilU] + Profiles.values(), 
                                  names=','.join(['U'] + Profiles.keys()))
    mlab.rec2csv(datatable, savefile, delimiter="\t")


def print_total_fluxes():
    hbeta = Fluxes["H__1__4861A"]
    print "H beta line flux: "
    print "Line fluxes relative to H beta = 100:"

    def sortkey(cloudy_line):
        """
        Key for sorting the lines by wavelength
        """
        cloudy_line = cloudy_line.replace('_', ' ')
        # strip off first 4 characters and then only keep the first word
        wav_s = cloudy_line[4:].split()[0]
        unit = wav_s[-1]
        wav = float(wav_s[:-1])
        if unit == 'A':
            wav *= 1.e-8
        elif unit == 'm':
            wav *= 1.e-4
        else:
            raise ValueError, "%s is neither A nor m" % (unit)
        return wav

    for emline in sorted(emlines, key=sortkey):
        flux = Fluxes[emline]
        em = emline.replace('_', ' ')
        print '|'.join(['', em[:-5], em[-5:], "%g" % (100.0*flux/hbeta), ''])



def calculate_profiles():
    """
    Returns 3 dicts (each keyed by line ID): 

    Fluxes : total flux of each lines
    Profiles : 1D velocity profile of each line
    DataCubes : 3D position-position-velocity cube of each line

    """
    Profiles = dict()
    Fluxes = dict()
    DataCubes = dict()

    if cmd_args.gauss:
        cubes = fastcube_gauss.calculate_cubes(R, Mu, V, T, Ems, awt, validmask, NU, NY, NX, NK, inc_degrees, xmin, ymin, umin, dx, dy, du)
    else:
        cubes = fastcube.calculate_cubes(R, Mu, V, Ems, validmask, NU, NY, NX, NK, inc_degrees, xmin, ymin, umin, dx, dy, du)

    print "Summing datacubes to find total fluxes...."
    for iline, emline in enumerate(emlines):
        # Put the physical units into the volume element
        DataCubes[emline] = cubes[iline,:,:,:]  * r0**3 
        # 1-D line profile as function of velocity
        Profiles[emline] =  DataCubes[emline].sum(axis=-1).sum(axis=-1)
        # Total line flux (integral of line profile)
        Fluxes[emline] = DataCubes[emline].sum()
        print Fluxes
    print ".... finished"


    return Fluxes, Profiles, DataCubes





if __name__ == "__main__": 

    cmd_args = parse_command_args()

    # Set some local variables from the command line arguments
    NU = cmd_args.nvel
    inc_degrees = cmd_args.inc
    r0 = cmd_args.r0
    Rmax = cmd_args.Rmax
    isRebinned = True

    xmin, ymin, xmax, ymax = cmd_args.window
    dy = dx = cmd_args.pixelsize
    # Note that xmin is the left border of the first cell
    # while xmax is the right border of the last cell
    NX = int((xmax - xmin)/dx)
    NY = int((ymax - ymin)/dy)
    print "Datacube size: (NU, NY, NX) = ({}, {}, {})".format(NU, NY, NX)


    if isRebinned:
        emlines, Mu, R, V, T, Ems, awt, validmask = read_rebinned_models()
    else:
        raise NotImplementedError

    #The model


    cosi = np.cos(np.radians(inc_degrees))
    sini = np.sin(np.radians(inc_degrees))

    # All arrays begin with capital letters
    phimin, phimax = 0.0, np.radians(360.0)
    thetamin, thetamax = 0.0, np.radians(90.0)

    Phi = np.linspace(phimin, phimax, num=cmd_args.nphi) #Lista de angulos en Phi
    NK = len(Phi)
    NJ = len(Mu)


    # Define the velocity bins
    umax = 50.0
    umin = -umax
    du = (umax - umin)/NU

    print "NI = {}, NJ = {}".format(len(R), len(Mu))
    print "Ems.shape = {}".format(Ems.shape)

    Fluxes, Profiles, DataCubes = calculate_profiles()

    PerfilU = np.linspace(umin, umax, NU)
    idstring = "nphi%(nphi)i-nvel%(nvel)i-inc%(inc)i" % (vars(cmd_args))
    savefile = "model_profile-%s.dat" % (idstring)


    write_lineprofile_table()

    print_total_fluxes()

    write_datacubes()

