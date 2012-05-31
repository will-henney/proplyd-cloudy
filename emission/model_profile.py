import numpy as np
import scipy
import sys, os, glob
import argparse
import pyfits

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
    Calculate 1-d line profile of Cloudy proplyd model. 
    If the --window option is used, 
    also calculate the 3-d position-position-velocity cube.
    """)

    parser.add_argument(
        "--modeldir", type=str, default=".",  
        help='Path to proplyd directory model')
    # parser.add_argument(
    #     "--ntheta", type=int, default=100,
    #     help='Number of angles theta')
    parser.add_argument(
        "--nphi", type=int, default=200,
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
        "--window", type=float, nargs=4, default=argparse.SUPPRESS,
        metavar=("X1", "Y1", "X2", "Y2"),
        help='If present, bounding box of view window (in units of r0) for the position-position-velocity cube')
    parser.add_argument(
        "--pixelsize", type=float, default=0.1, 
        help='Size of pixels (in units of r0) for the PPV cube')

    parser.add_argument(
        "--rebin", type=str, default=argparse.SUPPRESS,
        help='If present, a rebinned model will be loaded from "rebin-REBIN" instead of from the original "th??" directories')

    return parser.parse_args()


def read_original_models():
    """
    Read the emissivities, etc from a set od Cloudy models
    """

    # Find list of subdirectories with original model output
    thetadirs =  [
        os.path.basename(p) for p in 
        glob.glob(os.path.join(cmd_args.modeldir, "th??"))
        ]

    # For the moment, we use the existing thetas with no interpolation
    Theta = np.array([np.radians(float(s[2:])) for s in thetadirs])

    # The following arrays are defined for each angle:
    # Radius, Velocity, Emissivity
    Radius_vectors = []
    Velocity_vectors = []
    Emissivities_vectors = []
    # List of names of emission lines
    emlines = None
    for thetadir in thetadirs:
        # Read in the Cloudy model for each angle
        claudia.CloudyModel.indir = os.path.join(cmd_args.modeldir, thetadir)
        claudia.CloudyModel.outdir = claudia.CloudyModel.indir
        # find the last iteration
        modelinputpath = glob.glob(os.path.join(claudia.CloudyModel.indir, "*.in"))[-1]
        # extract only the basename of the file, without .in extension
        modelname = os.path.splitext(os.path.basename(modelinputpath))[0]
        print "Indir: ", claudia.CloudyModel.indir
        print "Modelname: ", modelname
        m = claudia.CloudyModel(modelname) 
        Radius = Rmax*r0 - m.ovr.depth # physical radius from proplyd center
        R = Radius/r0            # array of dimensionless radius
        i2 = len(R[R>=1.0])             # find index where R = 1

        # Find the list of emission lines if necesssary
        if emlines is None:
            emlines = m.em.dtype.names[1:] # miss off "depth"

        # Calculate physical variables as function of radius for this theta
        Density = 10**(m.ovr.hden)
        Emissivities = dict()
        for emline in emlines:
            Emissivities[emline] = 10**(m.em[emline])
        sound_speed = np.sqrt(3./5.)*m.pre.cadwind_kms # isothermal sound speed
        n0 = Density[i2]
        Velocity = sound_speed[i2] * n0 / (R**2 * Density) 

        # Add them to the lists of vectors 
        Radius_vectors.append(Radius)
        Velocity_vectors.append(Velocity)
        Emissivities_vectors.append(Emissivities)


    return emlines, Theta, Radius_vectors, Velocity_vectors, Emissivities_vectors

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
    Theta = np.arccos(read_fits("mu")[:,0])
    R = read_fits("R")[0,:]
    i2 = len(R[R<1.0])             # find index where R = 1
    Density2D = 10**read_fits("ovr-hden")
    sound_speed2D = np.sqrt(3./5.)*read_fits("pre-cadwind_kms")
    emfiles = [
        os.path.basename(p) for p in 
        glob.glob(
            os.path.join(
                cmd_args.modeldir,
                "rebin-%s" % (cmd_args.rebin),
                "em-*.fits"
                )
            )
        ]

    emlines = [ s.split("-")[1].split(".")[0] for s in emfiles ]
    Emissivities2D = dict([(emline, 10**read_fits("em-%s" % (emline))) for emline in emlines])
    Radius_vectors = []
    Velocity_vectors = []
    Emissivities_vectors = []
    for j, th in enumerate(Theta):
        n0 = Density2D[j, i2]
        Velocity = sound_speed2D[j, i2] * n0 / (R**2 * Density2D[j, :])
        Radius_vectors.append(R)
        Velocity_vectors.append(Velocity)
        Emissivities = dict(
            [(emline, Emissivities2D[emline][j,:]) for emline in emlines]
            )
        Emissivities_vectors.append(Emissivities)

    print  Velocity_vectors
    return emlines, Theta, Radius_vectors, Velocity_vectors, Emissivities_vectors




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
    for emline in emlines:
        # Line profile as function of velocity
        Profiles[emline] = np.zeros(NU)
        # Total line flux (integral of line profile)
        Fluxes[emline] = 0.0
        if hasDataCube:
            DataCubes[emline] = np.zeros((NU, NY, NX))


    for k in range(NK):
        cphi = np.cos(Phi[k])
        sphi = np.sin(Phi[k])
        kneg = max(0, k - 1)
        kpos = min(NK-1, k + 1)
        dphi = 0.5*(Phi[kpos] - Phi[kneg])
        for j in range(NJ):
            # for each angle
            Radius = Radius_vectors[j]
            Velocity = Velocity_vectors[j]
            Emissivities = Emissivities_vectors[j]

            ctheta = np.cos(Theta[j])
            stheta = np.sin(Theta[j])
            jneg = max(0, j - 1)
            jpos = min(NJ-1, j + 1)
            dmu = -0.5*(np.cos(Theta[jpos]) - np.cos(Theta[jneg]))
            NI = len(Radius)
            for i in range(NI):
                if Velocity[i] == -1.0:
                    continue    # Velocity = -1 indicates invalid radii 
                ineg = max(0, i - 1)
                ipos = min(NI-1, i + 1)
                dr = 0.5*(Radius[ipos] - Radius[ineg])
                dvol =  abs(dphi * dmu * (Radius[i]**2) * dr)
                u = -Velocity[i]*(sini*stheta*cphi + cosi*ctheta)
                iu = int(np.floor((u-umin)/du))
                if hasDataCube:
                    # Coordinates in plane of sky
                    x = (Radius[i]/r0)*(cosi*stheta*cphi + sini*ctheta)
                    y = (Radius[i]/r0)*stheta*sphi
                    ix = int(np.floor((x-xmin)/dx))
                    iy = int(np.floor((y-ymin)/dy))
                    # print "(iy, ix, iu) = ({}, {}, {})".format(iy, ix, iu)

                # assert iu >= 0 and iu < NU, "Index (%i) out of bounds [%i:%i] of velocity array. u = %.3g, umin, umax = %.3g, %.3g" % (iu, 0, NU-1, u, umin, umax)
                # for each emission lines
                for emline in emlines:
                    # add in to line profile
                    if iu >= 0 and iu < NU: 
                        Profiles[emline][iu] += dvol * Emissivities[emline][i]
                    # and add in to total flux
                    Fluxes[emline] += dvol * Emissivities[emline][i]
                    # ... and add in to datacube if required
                    if ( hasDataCube
                         and iu >= 0 and iu < NU 
                         and ix >= 0 and ix < NX 
                         and iy >= 0 and iy < NY ):
                        DataCubes[emline][iu, iy, ix] += dvol * Emissivities[emline][i]
    return Fluxes, Profiles, DataCubes





if __name__ == "__main__": 

    cmd_args = parse_command_args()

    # Set some local variables from the command line arguments
    NU = cmd_args.nvel
    inc_degrees = cmd_args.inc
    r0 = cmd_args.r0
    Rmax = cmd_args.Rmax
    hasDataCube = "window" in vars(cmd_args) # do we want a data cube or not?
    isRebinned = "rebin" in vars(cmd_args)   # are we using the rebinned data?

    if hasDataCube:
        xmin, ymin, xmax, ymax = cmd_args.window
        dy = dx = cmd_args.pixelsize
        # Note that xmin is the left border of the first cell
        # while xmax is the right border of the last cell
        NX = int((xmax - xmin)/dx)
        NY = int((ymax - ymin)/dy)
        print "Datacube size: (NU, NY, NX) = ({}, {}, {})".format(NU, NY, NX)


    if isRebinned:
        emlines, Theta, Radius_vectors, Velocity_vectors, Emissivities_vectors = read_rebinned_models()
    else:
        emlines, Theta, Radius_vectors, Velocity_vectors, Emissivities_vectors = read_original_models()


    #The model


    cosi = np.cos(np.radians(inc_degrees))
    sini = np.sin(np.radians(inc_degrees))

    # All arrays begin with capital letters
    phimin, phimax = 0.0, np.radians(360.0)
    thetamin, thetamax = 0.0, np.radians(90.0)

    Phi = np.linspace(phimin, phimax, num=cmd_args.nphi) #Lista de angulos en Phi
    print "Theta: ", Theta
    # Theta = np.linspace(thetamin, thetamax, num=cmd_args.ntheta) #Lista de angulos en Theta
    NK = len(Phi)
    NJ = len(Theta)


    # Define the velocity bins
    umax = max(abs(Velocity_vectors[0]))
    umin = -umax
    du = (umax - umin)/NU

    Fluxes, Profiles, DataCubes = calculate_profiles()

    PerfilU = np.linspace(umin, umax, NU)
    idstring = "nphi%(nphi)i-nvel%(nvel)i-inc%(inc)i" % (vars(cmd_args))
    savefile = "model_profile-%s.dat" % (idstring)


    write_lineprofile_table()

    print_total_fluxes()


    if hasDataCube:
        write_datacubes()

