
import numpy as np
import scipy.special as sp
#### SKEW-NORMAL DISTRIBUTION #############################
#
# Normal-distribution with a skewness parameter
#
#
# mu:       mean
# sigma :   standard deviation/width
# alpha:    skewness parameter
# amplitude: log(strength)
# c:        constant background
def skewed_normal(x, mu, sigma, alpha, amplitude, c=0.0):

    constants = np.exp(amplitude)/(np.sqrt(2.0*np.pi)*sigma**2.0)
    pdfnorm = np.exp(-(x-mu)**2.0/(2.0*sigma**2.0))
    cdfnorm = 1.0 + sp.erf(alpha*(x-mu)/(np.sqrt(2.0)*sigma))

    res = constants*pdfnorm*cdfnorm + c
    return res



### SKEW-NORMAL WITH LOG(SIGMA) and LOG(C)
#
# variation of the skew-normal distribution
# with logarithmic parameters for sigma and c
# such that I can implement a scale prior on both
#
#
def logskew(x, mu, sigma, alpha, amplitude, c=None):

    sigma = np.exp(sigma)
    if not c:
        c = 0.0
    else:
        c = np.exp(c)
    return skewed_normal(x, mu, sigma, alpha, amplitude, c)


def asymmetric_envelope(x, trise, amplitude, c=0.0):

    amplitude = np.exp(amplitude)

    tend = x[-1] + x[1]
    tstart = x[0] - x[1]

    trise = trise*(tend - tstart)

    tfall = (tend-tstart) - trise
    env = envelope(x, tstart, tend, trise, tfall, amplitude)
    env = env + c
    return env


def envelope(x, tstart, tend, trise, tfall, amplitude):

    rexp = -trise/(x - tstart)
    fexp = tfall/(x - tend)

    res = amplitude*np.exp(rexp + fexp)

    return res



def const(freq, a):
    return np.array([np.exp(a) for x in freq])

#### FAST-RISE EXPONENTIAL DECAY (FRED) PROFILE
#
# useful for light curves, not strictly for MLE fitting
#
# f(x) = a*exp(2*t1/t2)**(1/2) * exp(-t1/(time-ts) - (time-ts)/t2)
# a = normalization 
# ts = pulse start time
# t1 = characteristic of burst rise
# t2 = characteristic of burst decay
def fred(time, a, tau1, tau2, c = 0.0, t0 = 0.5):

    #print("a: " + str(a))
    #print("tau1: " + str(tau1))
    #print("tau2: " + str(tau2))
    #print("len time: " + str(len(time)))
    tdiff = time - (time[0]) - t0
    #print("len tdiff: " + str(len(tdiff)))
    dt = tdiff[1]-tdiff[0]
    ts = stepfun(time, time[0]+t0+dt, max(time)+dt)
    #tdiff = time - t0
    #print("tdiff: " + str(tdiff[0]))
    e1 = (-tau1/(tdiff)) - ((tdiff)/tau2)
    e2 = (np.exp(2.0*(tau1/tau2)**0.5))
    counts = a*np.exp(e1)*e*ts + c
    cmod = []
    for co in counts:
        if isnan(co) or co == np.inf:
           cmod.append(c)
        else:
           cmod.append(co)
    #print('funcval in FRED: ' + str(counts))

    return cmod

def combfred(x, norm1, norm2, norm3, norm4, tau11, tau12, tau13, tau14, tau21, tau22, tau23, tau24, t01, t02, t03, t04, c):

        fred1 = fred(x, norm1, tau11, tau21, 1.0, t01)
        fred2 = fred(x, norm2, tau12, tau22, 1.0, t02)
        fred3 = fred(x, norm3, tau13, tau23, 1.0, t03)
        fred4 = fred(x, norm4, tau14, tau24, 1.0, t04)
        fredsum = np.array(np.array(fred1) + np.array(fred2) + np.array(fred3) + np.array(fred4))# + c

        return fredsum

def stepfun(x, a, b):
    y = np.zeros(len(x))
    x = np.array(x)
    minind = x.searchsorted(a)
    maxind = x.searchsorted(b)
    y[minind:maxind] = 1.0
    return y

def gaussian(x, mean, scale, norm, c=0.0):
    norm = np.exp(norm)
    c = np.exp(c)
    g = norm*np.exp(-(x-mean)**2.0/(2.0*scale**2.0)) + c
    return g


