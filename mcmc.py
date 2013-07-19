
from __future__ import print_function
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.ticker import MaxNLocator
import cPickle as pickle
import copy

import numpy as np
import scipy.optimize
import scipy.stats
import math


import generaltools as gt
import posterior


class MarkovChainMonteCarlo(object):

    def __init__(self, x, y, lpost, logfile='test', datatype=None):

        ## dependent and independent variable
        self.x = x
        self.y = y

        ## datatype: if data is a periodogram or a GP, then I can use
        ## specialized methods:
        self.datatype = datatype

        ## posterior probability density function
        self.lpost = lpost

        ## write out messages and results in this logfile:
        self.logfile = gt.TwoPrint(test + '.dat')

    ### perform a maximum a posteriori fit on the data
    def map_fit(self, ain, fitmethod='bfgs'):

        if self.datatype in ['ps', 'power', 'powerspectrum', 'per', 'periodogram']:
            ps = 
            fitspec = mle.PerMaxLike(


