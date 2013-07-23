
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

    def __init__(self, x, y, lpost, logfile='test', namestr='test', datatype=None):

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

        ## string to use for filenames
        self.namestr = namestr
 
        return


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
    def _rhat(self):


        self.logfile("Computing Rhat. The closer to 1, the better!")

        rh = []

        ### loop over parameters ###
        for i,k in enumerate(self.popt):

            ### pick parameter out of array
            tpar = np.array([t[i] for t in self.flatchain])

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
            self.logfile("The Rhat value for parameter " + str(i) + " is: " + str(rh[i]) + ".")

            if rh[i] > 1.2:
                self.logfile("*** HIGH Rhat! Check results! ***")
            else:
                self.logfile("Good Rhat. Hoorah!")


        return rh


    def _quantiles(self):

        ### empty lists for quantiles
        ci0, ci1 = [], []

        ### loop over the parameters ###
        for i,k in enumerate(self.popt):

            print("I am on parameter: " + str(i))

            ### pick parameter out of array
            tpar = np.array([t[i] for t in self.flatchain])
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

    def check_convergence(self):


        ### compute Rhat for all parameters
        rh = self._rhat()
        self.rhat = rh

        plt.scatter(rh, np.arange(len(rh))+1.0 )
        plt.axis([0.1,2,0.5,0.5+len(rh)])
        plt.xlabel("R_hat")
        plt.ylabel("Parameter")
        plt.title('Rhat')
        plt.savefig(self.namestr + '_rhat.ps')
        plt.close()


        ### compute 80% quantiles
        ci0, ci1 = self._quantiles()


        ### set array with colours
        ### make sure there are enough colours available
        colours_basic = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
        cneeded = int(math.ceil(len(ci0[0])/7.0))
        colours = []
        for x in range(cneeded):
            colours.extend(colours_basic)

        ### plot 80% quantiles        if self.plot:
            plt.plot(0,0)            plt.axis([-2, 2, 0.5, 0.5+len(ci0)])
            for j in range(self.nchain):
                plt.hlines(y=[m+(j)/(4.0*self.nchain) for m in range(len(ci0))], xmin=[x[j] for x in ci0], xmax=[x[j] for x in ci1], colo
r=colours[j])
            #plt.hlines(y=[m+1.0+(1)/(4*self.nchain) for m in np.arange(len(ci0))], xmin=[x[1] for x in ci0], xmax=[x[1] for x in ci1], c
olor=colours[j])

            plt.xlabel("80% region (scaled)")
            plt.ylabel("Parameter")
            plt.title("80% quantiles")
            plt.savefig(self.namestr + "_quantiles.ps")
            plt.close()

    def mcmc_infer(self):

        if printobj:
            print = printobj
        else:
            from __builtin__ import print as print


        ### covariance of the parameters from simulations
        covsim = np.cov(self.flatchain)

        print("Covariance matrix (after simulations): \n")
        print(str(covsim))

        ### calculate for each parameter its (posterior) mean and equal tail
        ### 90% (credible) interval from the MCMC

        self.mean = map(lambda y: np.mean(y), self.flatchain)
        self.std = map(lambda y: np.std(y), self.flatchain)
        self.ci = map(lambda y: quantiles(y, prob=[0.05, 0.95]), self.flatchain)


        ### print to screen
        print("-- Posterior Summary of Parameters: \n")
        print("parameter \t mean \t\t sd \t\t 5% \t\t 95% \n")
        print("---------------------------------------------\n")
        for i in range(len(self.popt)):
            print("chain[" + str(i) + "] \t " + str(self.mean[i]) + "\t" + str(self.std[i]) + "\t" + str(self.ci[i][0]) + "\t" + str(self.ci[i][1]) + "\n" )

        #np.random.shuffle(self.mcall)

        ### produce matrix scatter plots

        ### number of parameters
        N = len(self.popt)
        print("N: " + str(N))
        n, bins, patches = [], [], []

        if self.plot:
            fig = plt.figure(figsize=(15,15))
            plt.subplots_adjust(top=0.925, bottom=0.025, left=0.025, right=0.975, wspace=0.2, hspace=0.2)
            for i in range(N):
                for j in range(N):
                    xmin, xmax = self.flatchain[j][:1000].min(), self.flatchain[j][:1000].max()
                    ymin, ymax = self.flatchain[i][:1000].min(), self.flatchain[i][:1000].max()
                    ax = fig.add_subplot(N,N,i*N+j+1)
                    #ax.axis([xmin, xmax, ymin, ymax])
                    ax.xaxis.set_major_locator(MaxNLocator(5))
                    ax.ticklabel_format(style="sci", scilimits=(-2,2))

                    #print('parameter ' + str(i) + ' : ' + str(self.popt[i]))

                    if i == j:
                        #pass
                        ntemp, binstemp, patchestemp = ax.hist(self.flatchain[i][:1000], 30, normed=True, histtype='stepfilled')
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
                        ax.scatter(self.flatchain[j][:1000], self.flatchain[i][:1000], s=7)
                        ### then add contours

#                        np.random.shuffle(self.mcall)

                        xmin, xmax = self.flatchain[j][:1000].min(), self.flatchain[j][:1000].max()
                        ymin, ymax = self.flatchain[i][:1000].min(), self.flatchain[i][:1000].max()


                        ### Perform Kernel density estimate on data
                        try:
                            X,Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
                            positions = np.vstack([X.ravel(), Y.ravel()])
                            values = np.vstack([self.flatchain[j][:1000], self.flatchain[i][:1000]])
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



            plt.savefig(self.namestr + "_scatter.png", format='png')
            plt.close()
        return




def MetropolisHastings(MarkovChainMonteCarlo, object):

    def __init__(self, x, y, lpost, logfile='test', datatype=None, emcee=True, pdist='mvn'):

        ## use emcee?
        self.emcee = emcee


        ## in theory, one could define a whole range of proposal distributions to be used
        ## at the moment, only a multivariate normal works, but should be easy to extend
        if pdist.lower() in ['mvn', 'multivariate normal', 'normal', 'gaussian']:
             self.pdist = np.random.mutlivariate_normal
        else:
             raise Exception('Proposal distribution not recognized!')

        MarkovChainMonteCarlo.__init__(x, y, lpost, logfile=logfile, datatype=datatype)
        return


    def run_mcmc(self, popt=None, cov=None, nchain=200, niter=100, burnin=100, a=2.0, namestr='test')


        ## if parameter burnin is smaller than one, interpret it as a fraction of niter
        if burnin < 1.0:
            self.burnin = burnin*niter
        ## else assume it's the actual number of samples to be thrown away
        else:
            self.burnin = burnin

        if not popt and not cov:
            try:
                self.popt
                self.cov
            except AttributeError:
                raise Exception('No parameter set and covariance matrix given. Aborting ...')

        else:
            self.popt = popt
            self.cov = cov

 
        p0 = [self.pdist(self.popt,self.pcov) for i in xrange(nchain)]

        allsamplers = []

        for nc in xrange(nchain):

           sampler = MHChain(self.popt, self.cov, self.niter, burnin=self.burnin, pdist=self.pdist, emcee=self.emcee)
           sampler.run_chain(t0=p0[nc])

           sampler.run_diagnostics(namestr=namestr+'_chain'+str(nc)+".dat")
           allsamplers.append(sampler)


        self.chain = [a.chain for a in allsamplers]
        flatchain()

    def flatchain(self):

        self.flatchain = list(self.chain[0])
        for c in self.chain[1:]:
            self.flatchain.extend(c)

        return










###############################
###############################
###############################



class MarkovChain(object):

    def __init__(self, popt, pcov, niter, burnin=0.5, pdist='mvn'): 

        if pdist.lower() in ['mvn', 'multivariate normal', 'normal', 'gaussian']:
             self.pdist = np.random.multivariate_normal

        self.popt = popt
        self.pcov = pcov

        self.niter = niter

        if burnin < 1.0:
            self.burnin = burnin*niter
            self.allsamples = niter
        else:
            self.burnin = burnin
            self.allsamples = burnin + niter


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

        for i,th in enumerate(self.chain[0]):
            #print "i: " + str(i)
            ts = np.array([t[i] for t in self.chain])
  
            #print "(i*3)+1: " + str((i*3)+1)
            ### plotting time series ###
            p1 = plt.subplot(len(self.popt), 3, (i*3)+1)
            p1 = plt.plot(ts)
            plt.axis([0, len(ts), min(ts), max(ts)])
            plt.xlabel("Number of draws")
            plt.ylabel("parameter value")
            plt.title("Time series for parameter " + str(paraname[i]) + ".")

            ### make a normal distribution to compare to
            #tsnorm = [np.random.normal(self.pcov[i], self.terr[i]) for x in range(len(ts)) ]

            p2 = plt.subplot(len(self.popt), 3, (i*3)+2)

            ### plotting histogram
            p2 = count, bins, ignored = plt.hist(ts, bins=10, normed=True)
            bnew = np.arange(bins[0], bins[-1], (bins[-1]-bins[0])/100.0)
            p2 = plt.plot(bnew, 1.0/(self.terr[i]*np.sqrt(2*np.pi))*np.exp(-(bnew - self.popt[i])**2.0/(2.0*self.terr[i]**2.0)), linewidth=2, color='r')
            plt.xlabel('value of ' + str(paraname[i]))
            plt.ylabel('probability')
            plt.title("Histogram for parameter " + str(paraname[i]) + ".")


            nlags = 30

            p3 = plt.subplot(len(self.popt), 3, (i*3)+3)
            acorr = gt.autocorr(ts,nlags=nlags, norm=True)
            p3 = plt.vlines(range(nlags), np.zeros(nlags), acorr, colors='black', linestyles='solid')
            plt.axis([0.0, nlags, 0.0, 1.0])
        #plt.show()
        plt.savefig(namestr  + "_diag.png", format='png',orientation='landscape')
        plt.close()
        return


class MHChain(MarkovChain, object):

    def __init__(self, popt, pcov, niter, burnin=0.5, pdist="mvn", emcee=True):
        self.emcee = emcee
        MarkovChain.__init__(popt, pcov, niter, burnin, pdist)


    def run_chain(self, t0=None, niter = None, burnin = None):

        if not niter == None:
            self.niter = niter
        if not burnin == None: 
            if burnin < 1.0:
                self.burnin = burnin*niter
                self.allsamples = self.niter
            else:
                self.burnin = burnin
                self.allsamples = burnin + self.niter



        if self.emcee:

            sampler = emcee.MHSampler(self.pcov, dim=ndim, lpostfn=lpost, args=[False])
            pos, prob, state = sampler.run_mcmc(p0[nc], self.burnin)
            sampler.reset()
 
            sampler.run_mcmc(pos, niter, rstate0=state)
 
            self.chain = sampler.chain
            self.lnprobability = sampler.lnprobability


        else:


            ### set up array
            ttemp, logp = [], []
            ttemp.append(self.t0)
            #lpost = posterior.PerPosterior(self.ps, self.func)
            logp.append(self.lpost(self.t0, neg=False))
     
            #print "self.popt: " + str(self.t0)
            #print "self.pcov: " + str(self.pcov)
     
            #print("np.arange(self.niter-1)+1" +  str(np.arange(self.niter-1)+1))
            for t in np.arange(self.allsamples-1)+1:
                print("cov: " + str(self.pcov))
     
                tprop = dist(ttemp[t-1], self.pcov)
                print("tprop: " + str(tprop))
     
                pprop = self.lpost(tprop)#, neg=False)
                #pprop = lpost(tprop, self.func, ps)
                print("pprop: " + str(pprop))
     
                #logr = logp[t-1] - pprop
                logr = pprop - logp[t-1]
                logr = min(logr, 0.0)
                r= np.exp(logr)
                update = choice([True, False], size=1, weights=[r, 1.0-r])
                print("update: " + str(update))
     
                if update:
                    ttemp.append(tprop)
                    logp.append(pprop)
                    if t > self.discard:
                        accept = accept + 1
                else:
                    ttemp.append(ttemp[t-1])
                    logp.append(logp[t-1])
     
            self.chain = ttemp[self.burnin+1:]
            self.lnprobability = logp[self.burnin+1:]
 
        return 



