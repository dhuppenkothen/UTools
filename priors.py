###### PRIORS ########
#
# Definitions for various prior functions
#
#
#
#
#
#
import scipy.stats
import numpy

def flatprior(x, xmin=None, xmax = None):

    if xmin is None and xmax is None:
        return gaussprior(x, 0.0, 1000.0)
    elif xmin is None and not xmax is None:
        xswitch = (x <= xmax)
        if xswitch:
            return gaussprior(x, 0.0, 1000.0)/2.0
        else:
            return 0.0

    elif xmax is None and not xmin is None:
        xswitch = (x >= xmin)
        if xswitch:
            return gaussprior(x, 0.0, 1000.0)/2.0
        else:
            return 0.0


    else:
        xswitch = (xmin <= x and x <= xmax)
        if xswitch:
            return 1.0/(xmax - xmin)
        else:
            return 0.0

def gaussprior(x, mean, sigma):
    px = scipy.stats.norm.pdf(x, mean, sigma)
    return px


