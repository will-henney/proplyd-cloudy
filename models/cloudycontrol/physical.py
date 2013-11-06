"""
Physical conditions in a Cloudy model

density, radius, etc

"""

from textwrap import dedent

def abundances(variant="Orion"):
    if variant == "Orion":
        abun_string = dedent("""\
                             * Orion nebula abundances plus Orion dust
                             abundances H II region no grains
                             grains Orion
                             grains PAH
                             set PAH "H" // Only have PAH in the neutral gas
                             """)
    elif variant == "FastOrion":
        abun_string = dedent("""\
                             * Orion nebula abundances plus simplified Orion dust
                             abundances H II region no grains
                             grains Orion single
                             """)
    elif variant == "TsamisLV2":
        abun_string = dedent("""\
                             * Tsamis abundances plus full Orion dust
                             abundances H II region no grains
                             element scale factor helium 1.09
                             element scale factor carbon 2.19
                             element scale factor nitrogen 1.02
                             element scale factor oxygen 2.48
                             element scale factor neon 3.16
                             element scale factor sulphur 0.68
                             element scale factor chlorine 2.29
                             element scale factor argon 1.29
                             element scale factor iron 0.03
                             grains Orion
                             grains PAH
                             set PAH "H" // Only have PAH in the neutral gas
                             """)
    elif variant == "Esteban":
        abun_string = dedent("""\
                             * Esteban et al 2004 (t2=0.002) abundances plus full Orion dust
                             abundances H II region no grains
                             element scale factor helium 1.02
                             element scale factor carbon 0.87
                             element scale factor nitrogen 0.76
                             element scale factor oxygen 1.12
                             element scale factor neon 1.86
                             element scale factor sulphur 1.66
                             element scale factor chlorine 2.88
                             element scale factor argon 1.38
                             element scale factor iron 0.33
                             grains Orion
                             grains PAH
                             set PAH "H" // Only have PAH in the neutral gas
                             """)
    elif variant == "HST1":
        abun_string = dedent("""\
                             * Further tweaked abundances to try and improve fit
                             abundances H II region no grains
                             element scale factor helium 1.0
                             element scale factor carbon 0.85
                             element scale factor nitrogen 1.0
                             element scale factor oxygen 0.5
                             element scale factor neon 0.6
                             element scale factor sulphur 0.95
                             element scale factor chlorine 2.0
                             element scale factor argon 0.6
                             element scale factor iron 0.15
                             grains Orion
                             grains PAH
                             set PAH "H" // Only have PAH in the neutral gas
                             """)
    elif variant == "HST10":
        abun_string = dedent("""\
                             * Further tweaked abundances to try and improve fit
                             abundances H II region no grains
                             element scale factor nitrogen 0.59
                             element scale factor oxygen 1.07
                             element scale factor sulphur 0.56
                             element scale factor neon 2.4
                             element scale factor iron 0.2
                             element scale factor argon 0.67
                             grains Orion
                             grains PAH
                             set PAH "H" // Only have PAH in the neutral gas
                             """)
    elif variant == "HST10nd":
        abun_string = dedent("""\
                             * Further tweaked abundances to try and improve fit
                             abundances H II region no grains
                             element scale factor nitrogen 0.59
                             element scale factor oxygen 1.07
                             element scale factor sulphur 0.56
                             element scale factor neon 2.4
                             element scale factor iron 0.2
                             element scale factor argon 0.67
                             """)
    elif variant == "Tweak01":
        abun_string = dedent("""\
                             * Further tweaked abundances to try and improve fit
                             abundances H II region no grains
                             element scale factor nitrogen 0.59
                             element scale factor oxygen 0.95
                             element scale factor sulphur 0.46
                             element scale factor neon 2.6
                             element scale factor iron 0.1
                             element scale factor argon 0.67
                             grains Orion
                             grains PAH
                             set PAH "H" // Only have PAH in the neutral gas
                             """)
    elif variant == "Tweak02":
        abun_string = dedent("""\
                             * Further tweaked abundances to try and improve fit
                             abundances H II region no grains
                             element scale factor nitrogen 0.59
                             element scale factor oxygen 0.95
                             element scale factor sulphur 0.46
                             element scale factor neon 3.0
                             element scale factor iron 0.08
                             element scale factor argon 0.67
                             grains Orion
                             grains PAH
                             set PAH "H" // Only have PAH in the neutral gas
                             """)
    elif variant == "Tweak03":
        abun_string = dedent("""\
                             * Further tweaked abundances to try and improve fit
                             abundances H II region no grains
                             element scale factor nitrogen 0.59
                             element scale factor oxygen 1.00
                             element scale factor sulphur 0.46
                             element scale factor neon 3.0
                             element scale factor iron 0.08
                             element scale factor argon 0.67
                             grains Orion
                             grains PAH
                             set PAH "H" // Only have PAH in the neutral gas
                             """)
    elif variant == "Tweak04":
        abun_string = dedent("""\
                             * Further tweaked abundances to try and improve fit
                             abundances H II region no grains
                             element scale factor nitrogen 0.5
                             element scale factor oxygen 0.9
                             element scale factor sulphur 0.4
                             element scale factor neon 3.8
                             element scale factor iron 0.08
                             element scale factor argon 0.67
                             grains Orion
                             grains PAH
                             set PAH "H" // Only have PAH in the neutral gas
                             """)
    elif variant == "Tweak05":
        abun_string = dedent("""\
                             """)
    elif variant == "Tweak06":
        abun_string = dedent("""\
                             """)
    elif variant == "Tweak07":
        abun_string = dedent("""\
                             """)
    elif variant == "Tweak08":
        abun_string = dedent("""\
                             """)
    elif variant == "Tweak09":
        abun_string = dedent("""\
                             """)
    elif variant == "Tweak10":
        abun_string = dedent("""\
                             """)
    else:
        raise NotImplementedError
    return abun_string
    
def proplyd(r0, den_Rmax, Rmax=9.0, W=30.0, x0=10.8, composition="Orion"):
    prop_string = "* Proplyd density law from customized dense_fabden.cpp\n"
    prop_string += "dlaw 53 %.3e %.3f %.4f %.1f %.1f\n" % (r0, Rmax, den_Rmax, W, x0)
    prop_string += abundances(composition)
    return prop_string

