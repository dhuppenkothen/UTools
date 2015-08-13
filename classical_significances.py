

import numpy as np
import time as tsys

def pavnosig(power, nspec, nsim=1.0e9, verbose=False):

    if power*nspec > 30000:
	if verbose:
            print("Probability of no signal too miniscule to calculate.")
        return 0.0

    else:
        fn = pavnosigfun(power, nspec, nsim, verbose)
        if verbose:
            print("Pr(Averaged power lies above P if no signal present) = %.4e" %fn)
        return fn



def pavnosigfun(power, nspec, nsim = 1.0e6, verbose=False):

    tst = tsys.clock()

    if nsim < 1.0e7:
        chisquare = np.random.chisquare(2*nspec, size=nsim)/nspec
        if verbose:
            print("The mean of the distribution is %f" %np.mean(chisquare))
            print("The variance of the distribution is %f" %np.var(chisquare))

        pval_ind = np.where(power < chisquare)[0]
        pval = len(pval_ind)/nsim
    else:
        pval = pavnosigfun_idl(power, nspec)

    tend = tsys.clock()
    if verbose:
        print("computation time %f" %(tend-tst))
    return pval


def pavnosigfun_idl(power, nspec):

    sum = 0.0
    m = nspec-1

    pn = power*nspec

    while m >= 0:

        #print("m %i" %m)
        s = 0.0
        for i in xrange(int(m)-1):
            #print("i %i" %i)
            #print("m-i: %f" %(m-i))
            s += np.log(float(m-i))
            #print("s: %f" %s)


        #print("s: %f" %s)
        logterm = m*np.log(pn/2.0) - pn/2.0 - s
        #print("logterm: %f" %logterm)
        term = np.exp(logterm)
        #print("term: %f" %term)
        ratio = sum/term
        #print("ratio: %f" %ratio)

        if ratio > 1.0e15:
            return sum

        sum += term
        m -= 1

    #print("The probability of observing P = %f under the assumption of noise is %.7f" %(power, sum))

    return sum
