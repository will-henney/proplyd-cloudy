"""
Miscellaneous cloudy commands
"""

def optimize(nomolecules=True, nolevel2=True):
    opt_string = "* Optimizations\n"
    if (nomolecules):
        opt_string += "no molecules\n"
    if (nolevel2):
        opt_string += "no level 2 lines\n"
    return opt_string

def stopping(efrac=0.05, temperature=5.0):
    stop_string = "*Stopping criteria\n"
    stop_string += "stop efrac %.3f\n" % (efrac)
    stop_string += "stop temperature %.3f linear\n" % (temperature)
    return stop_string

def iterate(iterations=None):
    if iterations is None:
        iter_string = "iterate\n"
    elif iterations == 0:
        iter_string = ""
    else:
        iter_string = "iterate %i\n" % (iterations)
    return iter_string
    
        
