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
    else:
        raise NotImplementedError
    return abun_string
    
def proplyd(r0, den_Rmax, Rmax=9.0, A=174.0, x0=10.8):
    prop_string = "* Proplyd density law from customized dense_fabden.cpp\n"
    prop_string += "dlaw 53 %.3e %.3f %.4f %.1f %.1f\n" % (r0, Rmax, den_Rmax, A, x0)
    prop_string += abundances("Orion")
    return prop_string

