## BAYESIAN ANALYSIS FOR PERIODOGRAMS
#
# 
#
#
# TO DO LIST:
# - add functionality for mixture models/QPOs to mlprior
# - add logging
# - add smoothing to periodograms to pick out narrow signals
# - add proposal distributions to emcee implementation beyond Gaussian
#

#!/usr/bin/env python

from __future__ import print_function
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.ticker import MaxNLocator
import cPickle as pickle
import copy
#import matplotlib
#matplotlib.use('png')


### GENERAL IMPORTS ###
import numpy as np
import scipy.optimize
from scipy.stats.mstats import mquantiles as quantiles
import scipy.stats
import time as tsys
import math


### New: added possibility to use emcee for MCMCs
try:
   import emcee
   import acor
   emcee_import = True
except ImportError:
   print("Emcee and Acor not installed. Using Metropolis-Hastings algorithm for Markov Chain Monte Carlo simulations.")
   emcee_import = False


### OWN SCRIPTS
import generaltools as gt
import lightcurve
import powerspectrum
import mle
import posterior


### Hack for Numpy Choice function ###
#
# Will be slow for large arrays.
#
#
# Input: - data= list to pick from
#        - weights = statistical weight of each element in data
#          if no weights are given, all choices are equally likely
#        - size = number of choices to generate (default: one)
#
# Output: - either single entry from data, chosen according to weights
#         - or a list of choices 
#
# Note that unlike numpy.random.choice, this function has no "replace"
# option! This means that elements picked from data will *always* be
# replaced, i.e. can be picked again!
#
#
def choice_hack(data, weights=None, size=None):



    #print "len(data): " + str(len(data))
    ### if no weights are given, all choices have equal probability
    if weights == None:
        weights = [1.0/float(len(data)) for x in range(len(data))]

    #print("weights: " + str(weights))
    #print "sum of Weights: " + str(sum(weights))
    if not np.sum(weights) == 1.0:
        if np.absolute(weights[0]) > 1.0e7 and sum(weights) == 0:
            weights = [1.0/float(len(data)) for x in range(len(data))]
        else:
            raise Exception("Weights entered do not add up to 1! This must not happen!")


    #print "Sum of weights: " + str(np.sum(weights))

    ### Compute edges of each bin
    edges = []
    etemp = 0.0
    for x,y in zip(data, weights):
       etemp = etemp + y
       edges.append(etemp)

    ### if the np.sum of all weights does not add up to 1, raise an Exception


    ### If no size given, just print one number
    if size == None:
        randno = np.random.rand()    

    ### Else make sure that size is an integer
    ### and make a list of random numbers
    try:
        randno = [np.random.rand() for x in np.arange(size)]
    except TypeError:
        raise TypeError("size should be an integer!")

    choice_index = np.array(edges).searchsorted(randno)
    choice_data = np.array(data)[choice_index]

    return choice_data



### See if cutting-edge numpy is installed so I can use choice
try:
    from numpy.random import choice
     ### if not, use hack
except ImportError:
    choice = choice_hack





### Compute log posterior ###
#
# This function computes the log- posterior
# probability density from the log-likelihood
# and the log-prior
#
# Returns: log-posterior probability density
#
def lpost(t0, func, ps):
    ### log-likelihood at parameter set t0
    mlogl = mle.maxlike(ps.freq, ps.ps, func, t0)
    ### log of prior distribution at t0
    priorval = mlprior(t0, func)
    ### log posterior is log-likelihood + log-prior
    mlpost = mlogl + priorval
    return mlpost


### Compute prior densities ###
#
# This function computes prior densities for
# all parameters and the whole parameter set
# 
# Returns the log-prior density
#
#
# NOTE: so far, this can only do power laws and bent
# power laws! At some point, I need to add functionality
# for mixture models and Lorentzians!
#
#
#
def mlprior(t0, func):

    ### allowed range for PL indices
    alim = [-1.0, 8.0]

    ### power law index always first value
    alpha = t0[0]
    ### palpha is True if it is within bounds and false otherwise
    ### then pr = pbeta*pgamma if palpha = True, 0 otherwise
    palpha = (alpha >= alim[0] and alpha <= alim[1])
    ### normalization always second parameter
    #beta = t0[1]
    pbeta = 1.0
    ### Poisson noise level always last parameter
    #gamma = t0[-1]
    pgamma = 1.0
    pr = palpha*pbeta*pgamma

    ### if we have a power law, we're done
    ### for a bent power law, there are two 
    ### more parameters:
    if func == mle.bpl:
    #    delta = t0[2]
        pdelta = 1.0
    #    eps = t0[3]
        peps = 1.0
        pr = pr*pdelta*peps

    #else:
    #   raise Exception("Function not defined. Will give all parameters a normal distribution!")

    ### compute logarithm of prior density
    if pr > 0:
        mlp = np.log(pr) 
    else:
        ### if prior density is zero, set it to a really small value
        ### to avoid logarithm errors
        mlp = -100.0 

    return mlp


##########################################
#########################################




###########################################################
###########################################################
###########################################################


###########################################################
#
# class MarkovChainMonteCarlo: make a sample of MCMCs
#
#
#
# TO DO: Add functionality for proposal distributions beyond a Gaussian
#
#
#
#
#
class MarkovChainMonteCarlo(object):

# ps: power spectrum (data or simulation)
# niter (int): number of iterations in each chain
# nchain (int): number of chains to run
# topt (list): best-fit set of parameters for MLE
# tcov (np.array): covariance matrix from MLE
# discard (int): fraction of chain to discard (standard: niter/2)
# func: function to be fitted (standard functions in mle module)
# paraname (optional): list of parameter names
# check_conv (bool): check for convergence?
# namestr (str): filename root for plots 
# use_emcee (bool): use emcee (True) or Metropolis-Hastings (False) for MCMCs?
#    def __init__(self, ps, topt, tcov, 
#                 covfactor=1.0, 
#                 niter=5000, 
#                 nchain=10, 
#                 discard=None, 
#                 func=mle.pl,  
#                 paraname = None, 
#                 check_conv = True, 
#                 namestr='test', 
#                 use_emcee=True,
#                 plot=True,
#                 printobj = None,
#                 m=1):


    def __init__(self, x, y, topt, tcov, 
                 covfactor=1.0,
                 niter=5000,
                 nchain=10,
                 discard=None,
                 lpost = mle.pl,
                 paraname = None,
                 check_conv = True,
                 namestr='test',
                 use_emcee=True,
                 plot=True,
                 printobj = None,
                 m=1):


        ### Make sure not to include zeroth frequency
#        self.ps = ps
        self.m = m

        self.x = x
        self.y = y

        if printobj:
            print = printobj
        else:
            from __builtin__ import print as print      

        self.plot = plot
        print("<--- self.ps len MCMC: " + str(len(self.x)))
        ### set of optimal parameters from MLE fitting
        self.topt = topt
        print("mcobs topt: " + str(self.topt))
        ### covariances of fitted parameters
        self.tcov = tcov*covfactor
        print("mcobs tcov: " + str(self.tcov))

        ### number of iterations for MCMC algorithm
        self.niter = niter
        ### number of MCMC chains to be computed
        self.nchain = nchain
        ### Error in the fitted parameters
        self.terr = np.sqrt(np.diag(tcov))
        ### function that was fitted
        self.lpost = lpost

        if discard == None:
            discard = math.floor(niter/2.0)


        mcall = []

        ### if emcee package is not installed, enforce Metropolis-Hastings implementation
        if emcee_import == False:
            print("Emcee not installed. Enforcing M-H algorithm!")
            use_emcee = False

        ### if emcee should be used, then use code below 
        ### loop over all chains
        for i in range(nchain):

            #t0 = topt + choice([2.0, 3.0, -3.0, -2.0], size=len(topt))*self.terr
		    
            ### set up MarkovChain object
            mcout = MarkovChain(niter = niter, topt = self.topt, tcov = self.tcov, lpost=self.lpost, paraname = paraname, emcee=use_emcee)
            ### create actual chain
            mcout.create_chain(self.x, self.y)
   
            ### make diagnostic plots
            mcout.run_diagnostics(namestr = namestr +"_c"+str(i), paraname=paraname)
 
            mcall.extend(mcout.theta)

        self.L = mcout.L
        mcall = np.array(mcall)

        ### check whether chains/walkers converged
        if check_conv == True:
            self.check_convergence(mcall, namestr, printobj = printobj)            

        ### transpose list of parameter sets so that I have lists with one parameter each
        self.mcall = mcall.transpose()

        ### make inferences from MCMC chain, plot to screen and save plots
        self.mcmc_infer(namestr=namestr, printobj = printobj)


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


    def check_convergence(self, mcall, namestr, printobj=None, use_emcee = True):


        if printobj:
            print = printobj
        else:
            from __builtin__ import print as print

        ### compute Rhat for all parameters
        rh = self._rhat(mcall, printobj)               
        self.rhat = rh        

        plt.scatter(rh, np.arange(len(rh))+1.0 )
        plt.axis([0.1,2,0.5,0.5+len(rh)])
        plt.xlabel("R_hat")
        plt.ylabel("Parameter")
        plt.title('Rhat')
        plt.savefig(namestr + '_rhat.ps')
        plt.close()


        ### compute 80% quantiles
        ci0, ci1 = self._quantiles(mcall)


        ### set array with colours
        ### make sure there are enough colours available
        colours_basic = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
        cneeded = int(math.ceil(len(ci0[0])/7.0))
        colours = []
        for x in range(cneeded):
            colours.extend(colours_basic)


        ### plot 80% quantiles
        if self.plot:
            plt.plot(0,0)
            plt.axis([-2, 2, 0.5, 0.5+len(ci0)])
            for j in range(self.nchain):
                plt.hlines(y=[m+(j)/(4.0*self.nchain) for m in range(len(ci0))], xmin=[x[j] for x in ci0], xmax=[x[j] for x in ci1], color=colours[j])
            #plt.hlines(y=[m+1.0+(1)/(4*self.nchain) for m in np.arange(len(ci0))], xmin=[x[1] for x in ci0], xmax=[x[1] for x in ci1], color=colours[j])

            plt.xlabel("80% region (scaled)")
            plt.ylabel("Parameter") 
            plt.title("80% quantiles")
            plt.savefig(namestr + "_quantiles.ps")
            plt.close()

    def mcmc_infer(self, namestr='test', printobj = None):

        if printobj:
            print = printobj
        else:
            from __builtin__ import print as print


        ### covariance of the parameters from simulations
        covsim = np.cov(self.mcall)

        print("Covariance matrix (after simulations): \n")
        print(str(covsim))

        ### calculate for each parameter its (posterior) mean and equal tail
        ### 90% (credible) interval from the MCMC

        self.mean = map(lambda y: np.mean(y), self.mcall)
        self.std = map(lambda y: np.std(y), self.mcall)
        self.ci = map(lambda y: quantiles(y, prob=[0.05, 0.95]), self.mcall)


        ### print to screen
        print("-- Posterior Summary of Parameters: \n")
        print("parameter \t mean \t\t sd \t\t 5% \t\t 95% \n")
        print("---------------------------------------------\n")
        for i in range(len(self.topt)):
            print("theta[" + str(i) + "] \t " + str(self.mean[i]) + "\t" + str(self.std[i]) + "\t" + str(self.ci[i][0]) + "\t" + str(self.ci[i][1]) + "\n" )

        #np.random.shuffle(self.mcall)

        ### produce matrix scatter plots

        ### number of parameters
        N = len(self.topt)
        print("N: " + str(N))
        n, bins, patches = [], [], []
 
        if self.plot:
            fig = plt.figure(figsize=(15,15))
            plt.subplots_adjust(top=0.925, bottom=0.025, left=0.025, right=0.975, wspace=0.2, hspace=0.2)
            for i in range(N):
                for j in range(N):
                    xmin, xmax = self.mcall[j][:1000].min(), self.mcall[j][:1000].max()
                    ymin, ymax = self.mcall[i][:1000].min(), self.mcall[i][:1000].max()
                    ax = fig.add_subplot(N,N,i*N+j+1)
                    #ax.axis([xmin, xmax, ymin, ymax])
                    ax.xaxis.set_major_locator(MaxNLocator(5))
                    ax.ticklabel_format(style="sci", scilimits=(-2,2))

                    #print('parameter ' + str(i) + ' : ' + str(self.topt[i]))
   
                    if i == j:
                        #pass
                        ntemp, binstemp, patchestemp = ax.hist(self.mcall[i][:1000], 30, normed=True, histtype='stepfilled')
                        n.append(ntemp)
                        bins.append(binstemp)
                        patches.append(patchestemp)
                        ax.axis([ymin, ymax, 0, max(ntemp)*1.2])
#                       ax.axis([xmin, xmax, 0, max(ntemp)*1.2])
   
                    else:
                        #ax = fig.add_subplot(N,N,i*N+j+1)

#                        ax.axis([xmin, xmax, ymin, ymax])

                        ax.axis([xmin, xmax, ymin, ymax])
                   #     np.random.shuffle(self.mcall)

                        ### make a scatter plot first
                        ax.scatter(self.mcall[j][:1000], self.mcall[i][:1000], s=7)
                        ### then add contours

#                        np.random.shuffle(self.mcall)

                        xmin, xmax = self.mcall[j][:1000].min(), self.mcall[j][:1000].max()
                        ymin, ymax = self.mcall[i][:1000].min(), self.mcall[i][:1000].max()
                        #print("xmin and xmax: " + str(xmin) + " and " + str(xmax))

                        ### Perform Kernel density estimate on data
                        try:
                            X,Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
                            positions = np.vstack([X.ravel(), Y.ravel()])
                            values = np.vstack([self.mcall[j][:1000], self.mcall[i][:1000]])
                            kernel = scipy.stats.gaussian_kde(values)
                            Z = np.reshape(kernel(positions).T, X.shape)
 
                            ax.contour(X,Y,Z,7)
                        except ValueError:
                            print("Not making contours.")

#        for i in range(N):
#            for j in range(N):
#                plt.subplot(N,N,i*N+j+1)
#                ax.xaxis.set_major_locator(MaxNLocator(5))
#                if i == j:
#                    plt.axis([min(bins[i]), max(bins[i]), 0.0, max(n[i]*1.2)])
#                    ax.xaxis.set_major_locator(MaxNLocator(5))
#                else:
#                    plt.axis([min(bins[j]), max(bins[j]), min(bins[i]), max(bins[i])])
#                    ax.xaxis.set_major_locator(MaxNLocator(4))
  


            plt.savefig(namestr + "_scatter.png", format='png')
            plt.close()
        return


#### POSTERIOR PREDICTIVE CHECKS ################
# 
# Note: fpeak is calculated in mle.PerMaxLike.compute_stats
# and can be found in dictionary self.pl_r or self.bpl_r
#
    ## nsim [int] = number of simulations
    ## dist [str] = distribution, one of
    ##      "exp": exponential distribution (=chi2_2), np.random.exponential
    ##       "chisquare": chi^2 distribution with df degrees of freedom
    ## df [int] = degrees of freedom for chi^2 distribution
    def simulate_periodogram(self, func=mle.pl, nsim=5000):

        ### number of simulations is either given by the user,
        ### or defined by the number of MCMCs run!
        nsim = min(nsim,len(self.mcall[0]))

 
        ### define distribution
#        if dist == "exp":
#             print("Distribution is exponential!")
#             noise = np.random.exponential(size = len(self.ps.freq))
#        elif dist == "chisquare":
#             noise = np.random.chisquare(2*df, size=len(self.ps.freq))/(2.0*df)
#        else:
#             raise Exception("Distribution not recognized")

        if self.m == 1:
            noise = np.random.exponential(size=len(self.x))
        else:
            noise = np.random.chisquare(2, size=len(self.x))/(2.0*self.m)

        ### shuffle MCMC parameters
        theta = np.transpose(self.mcall)
        #print "theta: " + str(len(theta))
        np.random.shuffle(theta)


        jump = int(np.floor(nsim/10))

        fper = []
        fps = []
        percount = 1.0

        perlist = [x*100.0 for x in range(10)]
        for x in range(nsim):

            if x in perlist:
                print(str(percount*10) + "% done ...")
                percount += 1.0
            ### extract parameter set
            ain = theta[x]
            ### compute model 'true' spectrum
            mpower = self.func(self.x, *ain)

            ### define distribution
#            if dist == "exp":
#                #print("Distribution is exponential!")
#                noise = np.random.exponential(size = len(self.ps.freq))
#            elif dist == "chisquare":
#                noise = np.random.chisquare(2*df, size=len(self.ps.freq))/(2.0*df)
#            else:
#                raiseException("Distribution not recognized")
            if self.m == 1:
                print("m = 1")
                noise = np.random.exponential(size=len(self.x))
            else:
                print("m = " + str(self.m))
                noise = np.random.chisquare(2*self.m, size=len(self.x))/(2.0*self.m)




            ### add random fluctuations
            mpower = mpower*noise

            ### save generated power spectrum in a PowerSpectrum object
            mps = powerspectrum.PowerSpectrum()
            mps.freq = self.x
            mps.ps = mpower
            mps.df = self.x[1] - self.x[0]
            mps.n = 2.0*len(self.x)
            mps.nphots = mpower[0]

            fps.append(mps)
        return fps




############################################
##### CODE BELOW IS A BAD IDEA! IGNORE THIS UNTIL MarkovChain CLASS DEFINITION
###########################################

#### MAKE A PARAMETER SET OBJECT ####
#
#
#
#
class ParameterSet(object):

    ### par = parameter vector (1D)
    ### parname = list with parameter names
    def __init__(par, parname = None):

        if parname == None:
            x = ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta', 'iota', 'lappa', 'lambda', 'mu']

        else:
            if not len(par) == len(parname):
                raise Exception("The length of the parameter array and the length of the list with parameter names doesn't match!")
            x = parname
 
        ### set attributes for parameter space
        for i,p in enumerate(par):
            setattr(self, x[i], p)


#### Don't have any methods for now ####


#########################################################
#########################################################
#########################################################

#### MAKE A MARKOV CHAIN OBJECT ###
#
# QUESTION: How can I make an object with variable
# parameters?
#
#
#
#  NEED TO THINK ABOUT HOW TO GET ATTRIBUTES!
#
class MarkovChain(object):


    def __init__(self, mcsuper = None, niter = 5000, topt = None, tcov =None, lpost = None, paraname=None, discard=None, emcee=True):

 
        self.niter = niter
        self.topt = topt
        self.tcov = tcov
        self.terr = np.sqrt(np.diag(tcov))
        self.t0 = topt + choice([2.0, 3.0, -3.0, -2.0], size=len(topt))*self.terr
        self.emcee = emcee

        self.lpost = lpost
        self.terr = np.sqrt(np.diag(tcov))
        if discard == None:
            self.discard = int(niter/2)
        else:
            self.discard = int(discard)
        if paraname == None:
            self.paraname = ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta', 'iota', 'lappa', 'lambda', 'mu']
        else:
            self.paraname = paraname

    ### set up MCMC chain
    ### possible distributions: 
    ###   - 'mvn': multi-variate normal (default)
    ###   - 'stt': student t-test
    def create_chain(self, x, y, topt=None, tcov = None, t0 = None, dist='mvn'):

        if not topt == None:
            self.topt = topt
        if not tcov == None:
            self.tcov = tcov
        if not t0 == None:
            self.t0 = t0


        ### set up distributions
        if dist=='mvn': 
             dist = np.random.multivariate_normal

        ### need to think about multivariate student-t distribution
        ### not available in numpy/scipy
        #elif dist == 'stt':
        #     dist = np.random.standard_t
        #     cov = cov/3.0

        ### set acceptance value to zero
        accept = 0.0

        if self.emcee:

            sampler = emcee.MHSampler(self.tcov, dim=ndim, lpostfn=lpost, args=[False])
            ### run burn-in phase and reset sampler
            pos, prob, state = sampler.run_mcmc(p0[i], niter)
            sampler.reset()

            ### run actual MCMCs
            sampler.run_mcmc(pos, niter, rstate0=state)

            ### list of all samples stored in flatchain
            self.theta = sampler.flatchain
            self.logp = sampler.lnprob



        else:

            ### set up array
            ttemp, logp = [], []
            ttemp.append(self.t0)
            #lpost = posterior.PerPosterior(self.ps, self.func)
            logp.append(self.lpost(self.t0, neg=False))       

            #print "self.topt: " + str(self.t0)
            #print "self.tcov: " + str(self.tcov)

            #print("np.arange(self.niter-1)+1" +  str(np.arange(self.niter-1)+1))
            for t in np.arange(self.niter-1)+1:
#               print("cov: " + str(self.tcov))

                tprop = dist(ttemp[t-1], self.tcov)
#               print("tprop: " + str(tprop))

                pprop = self.lpost(tprop)#, neg=False)
                #pprop = lpost(tprop, self.func, ps)
#               print("pprop: " + str(pprop))

                #logr = logp[t-1] - pprop
                logr = pprop - logp[t-1]
                logr = min(logr, 0.0)
                r= np.exp(logr)
                update = choice([True, False], size=1, weights=[r, 1.0-r])
#               print("update: " + str(update))

                if update:
                    ttemp.append(tprop)
                    logp.append(pprop)
                    if t > self.discard:
                        accept = accept + 1
                else:
                    ttemp.append(ttemp[t-1])
                    logp.append(logp[t-1])
        
            self.theta = ttemp[self.discard+1:]
            self.logp = logp[self.discard+1:]


                #logr = logp[t-1] - pprop
                #print("logr: " + str(logr))
                #r = logr #np.exp(logr)

                #samplechoice = choice([True, False], size=1, weights=[0.1, 0.9])

                #if r > 0.0:# or not samplechoice:
                #    ttemp.append(tprop)
                #    logp.append(pprop)
                #    if t > self.discard:
                #         accept = accept+1

                #else:
                #    print "r: " + str(r)
                #    newset = choice([[tprop, pprop], [ttemp[t-1], logp[t-1]]], size=1, weights=[np.exp(r), 1.0-np.exp(r)])
                #    ttemp.append(newset[0][0])
                #    logp.append(newset[0][1]) 
                #    if newset[0][0].all() == tprop.all() and t > self.discard:
                #        
                #        accept = accept+ 1
#
#                print "t: " + str(t)
#
#                ### draw value from proposal distribution
#                tprop = dist(ttemp[t-1], self.tcov)
#                #print "tprop: " + str(tprop)
#                ### calculate log posterior density
#                pprop = lpost(tprop, self.func, ps)
#                print "pprop: " + str(pprop)
#                print "old pprop: " + str(logp[t-1])
#                ### compute ratio of posteriors at new and old
#                ### locations in terms of log posteriors
#                logr = logp[t-1] - pprop
#                #logr = pprop - logp[t-1]
#                print "logr: " + str(logr)
#
#                logr = min(logr, 0.0)
#                r = np.exp(logr)
#                print "r: " + str(r)
#                if r > 1.01:
#                    print("tprop: " + str(tprop))
#                    print("pprop: " + str(pprop))
#                    print("log(r): " + str(logr))
#                    print("r: " + str(r))
#
#                    raise Exception("r > 1! This cannot be true!")
# 
#                ### with chance r, the number is updated 
#                update = choice([True, False], size=1, weights=[r,1.0-r])
#                #print "update" + str(update)
#
#
#                if update:
#                    ttemp.append(tprop)
#                    logp.append(pprop)
#                    if t > self.discard:
#                        accept = accept+1
#                else:
#                    ttemp.append(ttemp[t-1])
#                    logp.append(logp[t-1])
#      
            #print "discard: " + str(self.discard)
            #print "self.discard: " + str(self.discard)

        self.theta = ttemp[self.discard+1:]
        self.logp = logp[self.discard+1:]
        self.L = self.niter - self.discard
        #print "self.niter: " + str(self.niter)
        #print "self.L: " + str(self.L)
        self.accept = accept/self.L
        return

    def run_diagnostics(self, namestr=None, paraname=None, printobj = None):

        if printobj:
            print = printobj
        else:
            from __builtin__ import print as print

        print("Markov Chain acceptance rate: " + str(self.accept) +".")

        if namestr == None:
            print("No file name string given for printing. Setting to 'test' ...")
            namestr = 'test'

        if paraname == None:
           paraname = ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta', 'iota', 'lappa', 'lambda', 'mu']


        fig = plt.figure(figsize=(12,10))
        adj =plt.subplots_adjust(hspace=0.4, wspace=0.4)

        for i,th in enumerate(self.theta[0]):
            #print "i: " + str(i)
            ts = np.array([t[i] for t in self.theta]) 
  
            #print "(i*3)+1: " + str((i*3)+1)
            ### plotting time series ###
            p1 = plt.subplot(len(self.topt), 3, (i*3)+1)
            p1 = plt.plot(ts)
            plt.axis([0, len(ts), min(ts), max(ts)])
            plt.xlabel("Number of draws")
            plt.ylabel("parameter value")
            plt.title("Time series for parameter " + str(paraname[i]) + ".")

            ### make a normal distribution to compare to
            #tsnorm = [np.random.normal(self.tcov[i], self.terr[i]) for x in range(len(ts)) ]

            p2 = plt.subplot(len(self.topt), 3, (i*3)+2)


            ### plotting histogram
            p2 = count, bins, ignored = plt.hist(ts, bins=10, normed=True)
            bnew = np.arange(bins[0], bins[-1], (bins[-1]-bins[0])/100.0)
            p2 = plt.plot(bnew, 1.0/(self.terr[i]*np.sqrt(2*np.pi))*np.exp(-(bnew - self.topt[i])**2.0/(2.0*self.terr[i]**2.0)), linewidth=2, color='r')
            plt.xlabel('value of ' + str(paraname[i]))
            plt.ylabel('probability')
            plt.title("Histogram for parameter " + str(paraname[i]) + ".")


            nlags = 30

            p3 = plt.subplot(len(self.topt), 3, (i*3)+3)
            acorr = gt.autocorr(ts,nlags=nlags, norm=True)
            p3 = plt.vlines(range(nlags), np.zeros(nlags), acorr, colors='black', linestyles='solid')
            plt.axis([0.0, nlags, 0.0, 1.0])
        #plt.show()
        plt.savefig(namestr  + "_diag.png", format='png',orientation='landscape')
        plt.close()


##############################################################
    



