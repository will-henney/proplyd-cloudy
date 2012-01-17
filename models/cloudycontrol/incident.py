"""
The incident radiaton field for a Cloudy model

The following types of radiation are provided:

+ star(log_phiH, Tstar, log_g=4.0, log_Z=0.0, atmosphere='BB')

+ background(crboost=1.0)

+ nebula()

"""
from textwrap import dedent


_BLACKBODY_RADIATION = """\
* Photoionization equilibrium
black body, T=%(Tstar).0f K 
""" 
_WMBASIC_RADIATION = """\
* WMBasic (Pauldrach et al. 2001) non-LTE, line-blanketed, 
*  and wind-blanketed hot stars
table star wmbasic %(Tstar).0f %(log_g).3f %(log_Z).3f 
"""
_ATLAS_RADIATION = """\
* Atlas (Castelli & Kurucz 2004) LTE, plane-parallel, 
*  hydrostatic model atmospheres
table star atlas odfnew Z+%(log_Z).0f %(Tstar).3f %(log_g).3f
"""
_TLUSTY_RADIATION = """\
* Tlusty (Lanz & Hubeny 2003) non-LTE, line-blanketed, 
*  plane-parallel, hydrostatic O and B star SED
table star tlusty OBstar 3-dim %(Tstar).0f %(log_g).3f %(log_Z).3f
"""
_KURUCZ_RADIATION = """\
* Kurucz (1979) similar to Atlas but obsolete 
* - used only for comparison with Baldwin (1991)
table star kurucz %(Tstar).0f
"""

_atmosphere_patterns = dict(
    BB = _BLACKBODY_RADIATION,
    WM = _WMBASIC_RADIATION,
    AT = _ATLAS_RADIATION,
    TL = _TLUSTY_RADIATION,
    KU = _KURUCZ_RADIATION
    )

def star(log_phiH, Tstar, atmosphere='BB', log_g=4.0, log_Z=0.0):
    """
    A stellar spectrum with ionizing flux and effective temperature    
    
    Optionally, specify gravity, metallicity, and atmosphere model.

    Supported atmosphere models:
    "BB" : Blackbody [log_g, log_Z ignored]
    "WM" : WMBasic (Pauldrach et al. 2001) NLTE, wind
    "AT" : ATLAS (Castelli & Kurucz 2004) LTE, hydrostatic
    "TL" : Tlusty (Lanz & Hubeny 2003) NLTE, hydrostatic
    "KU" : Kurucz (1979) OBSOLETE [log_g, log_Z ignored]

    """
    assert atmosphere in _atmosphere_patterns.keys(), \
        "Atmosphere %s not recognised, must be one of %s." % \
        (atmosphere, ", ".join(_atmosphere_patterns.keys()))
    s = _atmosphere_patterns[atmosphere] % dict(Tstar=Tstar, log_g=log_g, log_Z=log_Z)
    s += "phi(H) %.2f\n" % (log_phiH)
    return s

def background(crboost=1.0):
    """
    Standard background radiation: diffuse interstellar + CMB + Cosmic Rays

    Optionally boost the cosmic rays by linear factor of crboost
    """
    return dedent("""\
                  * Background radiation that is always there
                  cosmic ray background %.4f linear
                  cmb
                  table ism // Diffuse interstellar field
                  """ % (crboost))

def nebula():
    """
    Diffuse emission from HII region

    TODO
    This will require taking the "save reflected continuum" output from another model
    """
    raise NotImplementedError
