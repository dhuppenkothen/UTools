import numpy as np
import mle
import scipy.misc


logmin = -100000
#### POSTERIOR DENSITY OBJECT
#
#
#
#
#
#
#
#
class Posterior(object):

    def __init__(self,x, y, func):
        self.x = x
        self.y = y

        ### func is a parametric model
        self.func = func

    def _const_prior(self, t0):
        pr = 1.0
        return pr



    ### simplest uninformative log-prior, define alternatives in subclasses!
    def logprior(self, t0):
        ### completely uninformative prior: flat distribution for all parameters
        pr = 1.0 

        return log(pr)

    ### use standard definition of the likelihood as the product of all 
    def loglikelihood(self, t0, neg=False):
        loglike = np.sum(np.array([self.func(x, t0) for x in self.x]))
        return loglike

#    def __call__(self, time, *t0):
    def __call__(self, t0, neg=False):
#        print("<-- neg: " + str(neg))
#        print("<--- t0: " + str(t0))
        lpost = self.loglikelihood(t0) + self.logprior(t0)
        #print("<--- logprior: " + str(self.logprior(t0)))

        #print("<--- loglike: " + str(self.loglikelihood(t0)))

        if neg == True:
            #print("<--- ### lpost: " + str(lpost))
            return lpost
        else:
            #print("<--- ### lpost: " + str(-lpost))
            return -lpost
    

class PerPosterior(Posterior):
    ### initialize posterior object with power spectrum 
    ### and broadband noise model function
    ### note: at the moment, only mle.pl and mle.bpl are supported
    def __init__(self, ps, func):
       self.ps=ps
       print('I am here!')
       Posterior.__init__(self,ps.freq[1:], ps.ps[1:], func)
       #super(Posterior,self).__init__(ps.freq[1:], ps.ps[1:], func)

    ### prior densities for PL model
    ### choose uninformative priors
    def pl_prior(self,t0):
        alim = [-8.0, 8.0]


        ### power law index always first value
        alpha = t0[0]
        ### palpha is True if it is within bounds and false otherwise
        ### then pr = pbeta*pgamma if palpha = True, 0 otherwise
        palpha = (alpha >= alim[0] and alpha <= alim[1])
        ### normalization always second parameter
        #beta = t0[1]
#        palpha = 1.0
        pbeta = 1.0

        pgamma = 1.0

        #if (alpha < alim[0]) or (alpha > alim[1]):
        #    return np.inf
        #else:
        return palpha*pbeta*pgamma

    ### prior densities for BPL model
    ### choose uninformative priors
    def bpl_prior(self,t0):
        pr0 = self.pl_prior([t0[2], t0[1], t0[-1]])

        delta =t0[3]
        pdelta = (delta >= min(self.ps.freq))

        peps = 1.0
        return pr0*pdelta*peps




    def plqpo_prior(self, t0):

        noise = t0[2]
        pnoise = scipy.stats.norm.pdf(noise, 0.6932, 0.2)/2.0
        if pnoise < 0.01:
            pnoise = 0.0

        gamma = t0[3]

        gamma_min = np.log(self.ps.df)
        gamma_max = np.log(nu0/2.0)

        pgamma = (gamma >= gamma_min and gamma <= gamma_max)

        norm = t0[4]
        pnorm = ( 20.0 >= norm >= -10.0)
     
        pr0 = 1.0

        try:
            nu = t0[5]
            pnu = scipy.stats.norm.pdf(nu, 93.0, 5)/0.08
            if pnu < 0.01:
                pnu = 0.0
            return pr0*pnoise*pgamma*pnorm*pnu
        except IndexError:
            return pr0*pnoise*pgamma*pnorm      

    def bplqpo_prior(self, t0):


        noise = t0[4]
        pnoise = scipy.stats.norm.pdf(noise, 0.6932, 0.2)/2.0
        if pnoise < 0.01:
            pnoise = 0.0

        gamma = t0[5]

        gamma_min = np.log(self.ps.df)
        gamma_max = np.log(nu0/2.0)

        pgamma = (gamma >= gamma_min and gamma <= gamma_max)

        norm = t0[6]
        pnorm = ( 50.0 >= norm >= -5.0)


        pr0 = 1.0
        return pr0*pgamma*pnorm


    def qpofixprior(self, t0):
        print("Using QPO fix prior")
        pr0 = self.pl_prior(t0[:3])
        #gamma = t0[3]
        #gamma_min = 3*self.ps.df
        #gamma_max = 50.0
        #pgamma = (gamma >= gamma_min and gamma <= gamma_max)
        pnorm = 1.0
        return pnorm*pr0


    def qpo_prior(self,t0):

        gamma = t0[0]
        nu0 = t0[2]
        norm = t0[1]

        gamma_min = np.log(self.ps.df)
        gamma_max = np.log(nu0/2.0)

        pgamma = (gamma >= gamma_min and gamma <= gamma_max)
        pnu0 = 1.0
        pnorm = 1.0

        return pgamma*pnu0*pnorm


    ### For now, assume combmod is broadband + QPO
    def combmod_prior(self, t0):
        if len(t0) in [5,6]:
            ### PL + QPO (+ noise)
            pr = self.pl_prior(t0[:3])*self.qpo_prior(t0[3:])
        elif len(t0) in [3,4]:
            pr = self._const_prior(t0[0])*self.qpo_prior(t0[1:])

        else:
            pr = self.bpl_prior(t0[:5])*self.qpo_prior(t0[5:])

        return pr
           

    def const_prior(self, t0):
        noise = t0[0]
        pmean = np.log(np.mean(self.ps.ps[1:]))
        #print("pmean: " + str(pmean))
        #print("pmean/10: " + str(pmean/5.0))
        #print("pmean*10.0: " + str(pmean*5.0))
        pnoise = scipy.stats.norm.pdf(noise, pmean, pmean/2.0)

#        pnoise = (noise <= pmean/5.0 and noise >= pmean*5.0)
        #print("pnoise: " + str(pnoise))
        return pnoise




    ### log of the prior
    ### actually, this is -log(prior)
    ### useful so that we can compute -log_posterior
    def logprior(self, t0):

        if self.func == mle.pl:
           mlp = self.pl_prior(t0)
        elif self.func == mle.bpl:
           mlp = self.bpl_prior(t0)
        elif self.func == mle.bpl2:
           mlp = self.bpl_prior(t0)
        elif self.func == mle.qpo or self.func.func_name == "lorentz":
           #print("Using QPO prior")
           mlp = self.qpo_prior(t0)
#        elif self.func == mle.const:
#           mlp = self._const_prior(t0)
        elif self.func.func_name == "combmod":
           mlp = self.combmod_prior(t0)
        elif self.func.func_name == "qpofix":
           mlp = self.qpofixprior(t0)
        elif self.func.func_name == "plqpo":
           mlp = self.plqpo_prior(t0)
        elif self.func.func_name == "bplqpo":
           mlp = self.bplqpo_prior(t0)
        elif self.func == mle.const:
           #print("Running constant prior")
           mlp = self.const_prior(t0)

        else:
           print("not running constant prior")
           mlp = 1.0

        if mlp > 0:
           return -np.log(mlp)
        else:
           return -logmin
           #return -logmin

    ### LOG - LIKELIHOOD
    ### actually, this is -log-likelihood (!)
    ### minimizing -logL is the same as maximizing logL
    ### so it all works out okay
    def loglikelihood(self,t0, neg=False):
        #print("self.ps.freq: " + str(self.ps.freq))
        #print("t0: " + str(t0))
        funcval = self.func(self.ps.freq, *t0)

        res = np.sum(np.log(funcval))+ np.sum(self.ps.ps/funcval)
#        print('res: ' + str(res))
#        print("type res: " + str(type(res)))
        if np.isnan(res):
            print("res is nan")
            res = -logmin
        elif res == np.inf or np.isfinite(res) == False:
            print("res is infinite!")
            res = -logmin


        return res

class LightcurvePosterior(Posterior):
    ### initialize posterior object with power spectrum 
    ### and broadband noise model function
    ### note: at the moment, only mle.pl and mle.bpl are supported
    def __init__(self, lc, func):
       self.lc=lc
       Posterior.__init__(self,lc.time, lc.countrate, func)

    ### at the moment, fit using only data, no prior
    def logprior(self, t0):
        lpr = -np.log(self._const_prior(t0))
        print("<-- log prior: " + str(lpr))
        return -np.log(self._const_prior(t0))

    def loglikelihood(self,t0, neg=False):
        print('t0: ' + str(t0))
        funcval_temp = self.func(self.lc.time, *t0)
    
        funcval = np.zeros(len(funcval_temp))  
        for i,f in enumerate(funcval_temp):
            if f == 0:
                funcval[i] = 1.0e-100
            #elif f == np.inf:
            #    funcval[i] = 1.0e100
            else:
                funcval[i] = f       
 
        print("<--- funcval: " + str(np.log(funcval)))
        #res = np.sum(funcval - self.lc.counts + self.lc.counts*(np.log(self.lc.counts) - np.log(funcval)))
        res = 2.0*np.sum(funcval - self.lc.counts*np.log(funcval))
        #res = sum((self.lc.counts - funcval)**2.0)
        print('res: ' + str(res))
        #print("type res: " + str(type(res)))
        if np.isnan(res):
            res = 0.0
        elif res == np.inf:
#            print("res is infinite!")
            res = -logmin

        return res


class StackPerPosterior(PerPosterior, object):

    def __init__(self, ps, func, m):
        self.m = m
        PerPosterior.__init__(self, ps, func)


    def loglikelihood(self,t0, neg=False):

        funcval = self.func(self.ps.freq, *t0)

#        res = np.sum(np.log(funcval))+ np.sum(self.ps.ps/funcval)

        res = 2.0*self.m*(np.sum(np.log(funcval))+ np.sum(self.ps.ps/funcval) + np.sum((2.0/float(2*self.m) - 1.0)*np.log(self.ps.ps)))
#        print('res: ' + str(res))
#        print("type res: " + str(type(res)))
        if np.isnan(res):
            #print("res is nan")
            res = -logmin
        elif res == np.inf:
            #print("res is infinite!")
            res = -logmin

        return res

### POSTERIOR FOR GAUSSIAN PROCESS REGRESSION #####
#
# This subclass defines a Posterior for a Gaussian Process
# with covariance function 'covariance', which has 'ncovar'
# parameters. 
# It allows for the additional specification of a function 
# 'fun' to be multiplied with the exponent of the red noise
# process. 
#
# NOTE: At the moment, log(prior) = 1
#
# created 06/06/13
#
#
class GaussProcPosterior(Posterior):

    def __init__(self,x,y, covariance, ncovar=None, fun=None):
        ### fun is the envelope function; if None, then fun=1.0
        if fun == None:
            self.fun = 1.0
        else:
            self.fun = fun

        self.ncovar = ncovar 
  
        ### save covariance function
        self.covar = covariance

        Posterior.__init__(self,x, y, covariance)


    def logprior(self, t0):
        gamma = t0[0]
        norm = t0[1]
        pgamma = (1.0 <= gamma and gamma <= 5.0)
        pnorm = 1.0
        mlp = pgamma*pnorm

        if mlp > 0:
           return -np.log(mlp)
        else:
           return -logmin


        return pgamma*pnorm


    def loglikelihood(self, t0):

        if not self.ncovar == None:
            theta_cov = t0[:self.ncovar]
            theta_fun = t0[self.ncovar:]
        else:
            theta_cov = t0

        if not self.fun == 1.0:
            logfun = np.log(self.fun(self.x, *theta_fun))
            diffy = np.log(self.y) - logfun
        else:
            diffy = self.y

        covar_matrix = self.covar(self.x, *theta_cov)
#        print('covar_matrix: ' + str(covar_matrix))       

#        print('t0: ' + str(theta_cov))     
        #chol = np.linalg.cholesky(covar_matrix)
        #w = np.diagonal(covar_matrix)
        #wlog = np.log(w)
        #logdetc = 2.0*np.sum(wlog)
 
        ldetc = np.linalg.slogdet(covar_matrix)
        logdetc = ldetc[0]*ldetc[1]

        invc = np.linalg.inv(covar_matrix)
        k = float(len(self.y))

        ## NOTE: this is the NEGATIVE log likelihood
        res =  0.5*k*np.log(2.0*np.pi) + 0.5*np.dot(np.dot(np.transpose(diffy), invc), diffy) + 0.5*logdetc
#        print('loglike1: ' + str(-0.5*np.dot(np.dot(diffy, invc), diffy)))
#        print('loglike2: ' + str(0.5*logdetc))


        if np.isnan(res):
            print("res is nan")
            res = -logmin
        elif res == np.inf or np.isfinite(res) == False:
            print("res is infinite!")
            res = -logmin


        return res




 
