import numpy as np
from scipy import fft
import copy 


def convolution(xdata, ydata, nlags=100, norm=True):

    xnorm = np.array(xdata) - np.mean(xdata)
    ynorm = np.array(ydata) - np.mean(ydata)

    xstd = np.std(xnorm)
    ystd = np.std(ynorm)

    lx = len(xnorm)
    ly = len(ynorm)

#    if not lx == ly:
#        raise Exception("Length of time series not equal!")

    xext = np.zeros(3*len(xnorm))
    yext = np.zeros(3*len(ynorm))

    xext[lx:2*lx] = xnorm
    #yext[ly:2*ly] = ynorm

    cor = []

    for i in range(2*lx):

        yext[(2*lx-i):(3*lx-i)] = ynorm

        corval = np.sum(xext*yext)/(xstd*ystd)
        cor.append(corval)

    if not nlags is None:
        cor = cor[lx-nlags:lx+nlags]

    if norm:
        cor = cor/max(cor)


    return cor


def crosscorrelation(xdata, ydata, nlags=100, norm=True):

    xnorm = np.array(xdata) - np.mean(xdata)
    ynorm = np.array(ydata) - np.mean(ydata)

    xstd = np.std(xnorm)
    ystd = np.std(ynorm)

    lx = len(xnorm)
    ly = len(ynorm)

#    if not lx == ly:
#        raise Exception("Length of time series not equal!")

    xext = np.zeros(3*len(xnorm))
    yext = np.zeros(3*len(ynorm))
    
    xext[lx:2*lx] = xnorm
    #yext[ly:2*ly] = ynorm

    cor = []

    for i in range(2*lx):

        yext[i:i+lx] = ynorm

        corval = np.sum(xext*yext)/(xstd*ystd)
        cor.append(corval)

    if not nlags is None:
        cor = cor[lx-nlags:lx+nlags]

    if norm:
        cor = cor/max(cor)


    return cor
 

def crossspectrum(xdata, ydata):

    xnorm = np.array(xdata) - np.mean(xdata)
    ynorm = np.array(ydata) - np.mean(ydata)

    xstd = np.std(xnorm)
    ystd = np.std(ynorm)

    lx = len(xnorm)
    ly = len(ynorm)

    if not lx == ly:
        raise Exception("Length of time series not equal!")

    xext = np.zeros(3*len(xnorm))
    yext = np.zeros(3*len(ynorm))

    xext[lx:2*lx] = xnorm
    #yext[ly:2*ly] = ynorm

    cor = []

    for i in range(2*lx):

        yext[i:i+lx] = ynorm

        corval = np.sum(xext*yext)
        cor.append(corval)

    cross = fft(cor)
    if norm:
        cor = cor/max(cor) 

    return cor
