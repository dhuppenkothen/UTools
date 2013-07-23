
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


    def run_mcmc(self, popt=None, cov=None, nchain=200, niter=100, burnin=100, a=2.0)


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

 
        p0 = [self.pdist(self.popt,self.tcov) for i in xrange(nchain)]

        allsamplers = []

        for nc in xrange(nchain):
            if self.emcee:

                 sampler = emcee.MHSampler(self.tcov, dim=ndim, lpostfn=lpost, args=[False])
                 pos, prob, state = sampler.run_mcmc(p0[nc], self.burnin)
                 sampler.reset()

                 sampler.run_mcmc(pos, niter, rstate0=state)

                 allsamplers.append(sampler)

            else:



## need to define MarkovChain object here, so I can run several of them :)









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

    def run_chain(self, niter = None, burnin = None):

        if not niter == None:
            self.niter = niter
        if not burnin == None: 
            if burnin < 1.0:
                self.burnin = burnin*niter
                self.allsamples = self.niter
            else:
                self.burnin = burnin
                self.allsamples = burnin + self.niter


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





