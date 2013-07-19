
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
    def map_fit(self, ain, bounds=None, fitmethod='bfgs', obs=True, **kwargs):

        if self.datatype in ['ps', 'power', 'powerspectrum', 'per', 'periodogram']:
            ps = powerspectrum.PowerSpectrum()
            ps.freq = self.x
            ps.ps = self.y
            ps.n = 2.0*len(ps.freq)
            ps.nphots = ps.ps[0]
            ps.df = ps.freq[1] - ps.freq[0] 
            fitspec = mle.PerMaxLike(ps, fitmethod=fitmethod, obs=obs)
            self.fitparams = fitspec.mlest(self.lpost, ain, bounds=bounts)

        elif self.datatype in ['gauss', 'gp', 'gaussproc', 'gaussian process']:
            fitspec = mle.GaussMaxLike(self.x, self.y, obs=obs, fitmethod=fitmethod)
            self.fitparams = fitspec.mlest(self.lpost, ain, func=kwargs['func'], obs=obs, bounds=bounds)

        else:
            fitspec = mle.MaxLikelihood(self.x, self.y, obs=obs, fitmethod=fitmethod)
            self.fitparams = fitspec.mlest(self.lpost, ain, bounds=bounds, obs=obs)

        self.popt = fitparams['popt'] 
        self.cov = fitparams['cov']

        return fitparams


    ### auxiliary function used in check_convergence
    ### computes R_hat, which compares the variance inside chains to the variances between chains
    def _rhat(self, mcall, printobj = None):


        if printobj:
            print = printobj
        else:
            from __builtin__ import print as print

        print("Computing Rhat. The closer to 1, the better!")

        rh = []

        ### loop over parameters ###
        for i,k in enumerate(self.topt):

            ### pick parameter out of array
            tpar = np.array([t[i] for t in mcall])

            ### reshape back into array of niter*nchain dimensions
            tpar = np.reshape(tpar, (self.nchain, len(tpar)/self.nchain))

            ### compute mean of variance of each chain

            #### THIS DOESN'T WORK FOR SOME REASON! TAKES VARIANCE OF EACH ELEMENT!!!
            ### CHECK THIS!
            sj = map(lambda y: np.var(y), tpar)
            W = np.mean(sj)

            ### compute variance of means of each chain
            mj = map(lambda y: np.mean(y), tpar)
            ### note: this assumes the discards
            B = np.var(mj)*self.L

            ## now compute marginal posterior variance
            mpv = ((float(self.L)-1.0)/float(self.L))*W + B/float(self.L)

            ### compute Rhat
            rh.append(np.sqrt(mpv/W))

            ### print convergence message on screen:
            print("The Rhat value for parameter " + str(i) + " is: " + str(rh[i]) + ".")

            if rh[i] > 1.2:
                print("*** HIGH Rhat! Check results! ***")
            else:
                print("Good Rhat. Hoorah!")


        return rh


    def _quantiles(self, mcall):

        ### empty lists for quantiles
        ci0, ci1 = [], []

        ### loop over the parameters ###
        for i,k in enumerate(self.topt):

            print("I am on parameter: " + str(i))

            ### pick parameter out of array
            tpar = np.array([t[i] for t in mcall])
            ### reshape back into array of niter*nchain dimensions
            tpar = np.reshape(tpar, (self.nchain, len(tpar)/self.nchain))

            ### compute mean of variance of each chain
            intv = map(lambda y: quantiles(y, prob=[0.1, 0.9]), tpar)

            ### quantiles will return a list with two elements for each
            ### chain: the 0.1 and 0.9 quantiles
            ### need to pick out these for each chain
            c0 = np.array([x[0] for x in intv])
            c1 = np.array([x[1] for x in intv])

            ### now compute the scale
            scale = np.mean(c1-c0)/2.0

            ### compute means of each chain
            mt = map(lambda y: np.mean(y), tpar)
            ### mean of means of all chains
            offset = np.mean(mt)

            ### rescale quantiles (WHY??)
            ci0.append((c0 - offset)/scale)
            ci1.append((c1 - offset)/scale)

        return ci0, ci1



