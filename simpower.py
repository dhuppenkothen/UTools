

import numpy as np
import scipy
import math
import argparse

import lightcurve
import powerspectrum


#### MULTIPLY LIGHT CURVES TOGETHER ##################
#
# Little functions that multiplies light curves of different 
# processes together.
# Base_lc can be any LightCurve object, even one of the three options given below
# Base_lc should be normalized to the desired mean, rednoise and qpo to 1 and the envelope
# to 1/mean(flux)
#
# base_lc [LightCurve] = base light curve to use --> should be the longest one
# env [list]= burst envelope, deterministic function
# rednoise [list] = red noise profile
# QPO [list] = quasi-periodic oscillation
#
# !!! IMPORTANT!!! MAKE SURE envelope, rednoise AND qpo ARE LISTS, NOT NUMPY ARRAYS!!!
#
#
#
#
def multiply_lightcurves(base_lc, envelope=None, rednoise=None, qpo=None):

    if envelope:
        base_lc.counts[:len(envelope)] = base_lc.counts[:len(envelope)]*envelope
    
    if rednoise:
        base_lc.counts = base_lc.counts*rednoise

    if qpo:
        base_lc.counts = base_lc.counts*qpo
   
    return


#############################################
#
# python implementation of Timmer+Koenig 1995
# simulations of red noise
#
#
#
class TimmerPS(powerspectrum.PowerSpectrum):

    ### fnyquist = number of points in light curve / 2.0
    ### dnu = frequency resolution, dnu = 1.0/length of time interval
    ### rms = fractional rms amplitude of the light curve
    ### nu0 = centroid frequency of Lorentzian (for QPO)
    ### gamma = width of Lorentzian (for QPO)
    ### pind = power law index
    ### pnorm = normalization of power law
    ### nphot = number of photons in light curve
    ### psshape = shape of power spectrum
    ### psmanual = put in array with chosen shape
    def __init__(self, fnyquist=4096.0, dnu=1.0, rms=1.0, nu0=None, gamma=None, pind=None, pnorm = None, nphot=None, psshape='plaw', psmanual=None):

	### make an empty PowerSpectrum object
	powerspectrum.PowerSpectrum.__init__(self,lc=None, counts=None)
	#### CREATE ARTIFICIAL POWER SPECTRUM
  
	### number of elements in lc/ps
	N = np.ceil(2.0*fnyquist/dnu)
        
        #print "N: " + str(N)
	### frequency array        
        self.n = N
	self.freq = np.arange(math.floor(N/2.0))*dnu + dnu 
	self.fnyquist = fnyquist
	self.dnu = dnu
        self.nphot = nphot


        ### turn rms into a variance of the log-normal light curve
        lnvar = np.log(rms**2.0 + 1.0)
        #print("Variance of log-normal light curve: " + str(lnvar))

	### make a shape for the power spectrum, depending on
	### psshape specified

	if psshape.lower() in ['flat', 'constant', 'white', 'white noise']:
            ## assume white noise power spectrum, <P> = N*sigma_ln
            s = np.array([self.n*lnvar for x in self.freq])          

	elif psshape.lower() in ['powerlaw', 'plaw']:
	    s = self.n*lnvar*(1.0/self.freq)**pind

### Don't do other shapes for now, until I need them
### CAREFUL: normalization of these is not right yet!
        elif psshape.lower() in ['qpo', 'lorentzian', 'periodic']:
            #print('I am here!')
            alpha = (gamma/math.pi)*dnu*N/2.0
            sold = alpha/((self.freq-nu0)**2.0 + gamma**2.0)
            snew = sold/sum(sold)
            #print('sum snew: ' + str(sum(snew)))
            s = (sold/sum(sold))*lnvar*fnyquist*self.n/self.dnu
#        elif psshape.lower() in ['w+p', 'combined plaw']:
#            s = np.array([rms**2.0+pnorm*(1/x)**2.0 for x in self.freq])

#        elif psshape.lower() in ['w+q', 'combined qpo']:
#            alpha = (sigma**2.0)*(gamma/math.pi)*dnu*N/2.0
#            s = 2.0 + nphot*alpha/((self.freq-nu0)**2.0 + gamma**2.0)

        elif psshape.lower() in ['manual', 'psmanual']:
        
            if not psmanual == None:
                #print(sum(psmanual/sum(psmanual))) 
                ### for now, assume variance normalization
                #s = (psmanual/sum(psmanual))*lnvar*fnyquist*self.n**2.0/2.0
                s = (psmanual/sum(psmanual))*lnvar*fnyquist*self.n**2.0/(self.dnu)

                #s = (psmanual/sum(psmanual))*self.n*(self.n/2.0)*lnvar

            else:
                raise Exception("No shape given!")

        #sfinal = np.insert(s, 0, 0)

        #print "len(s) : " + str(len(s))
        #print "type(s): " + str(type(s))
        ### first element is zero, that will be the number of photons
        #sfinal = np.insert(s, 0, 0.0)

        ### s is the power spectrum
        self.s = s #sfinal


    def makeFourierCoeffs(self):
        nphot = self.nphot
        N = self.n 
 
        a = np.zeros(len(self.s))
        x1 = np.random.normal(size=len(self.s))*(self.s/2.0)**0.5
        y1 = np.random.normal(size=len(self.s))*(self.s/2.0)**0.5


        ### S(fnyquist) is real
        y1[-1] = 0.0
 
        self.x1 = x1
        self.y1 = y1

        ### now make complex Fourier pair
        Fpos = [complex(re,im) for re, im in zip(x1,y1)]
        Fneg = [complex(re,-im) for re, im in zip(x1,y1)]

        #print "Fpos: " + str(Fpos[:5])
        #print "Fpos: " + str(Fpos[-5:])
        #print "Fneg: " + str(Fneg[:5])
        #print "Fneg: " + str(Fneg[-5:])

        Fpos.insert(0, (0+0j))


        Fneg.reverse()

        #print "Fneg: " + str(Fneg[:5])
        #print "Fneg: " + str(len(Fneg))
        #print "Fneg: " + str(Fneg[-5:])


        ### remove duplicate nyquist frequency and the nu=0 term
        #Fneg = Fneg[1:1+int(np.round((N-1)/2))]
        #Fneg = Fneg[:-1]
        #print "Fneg: " + str(len(Fneg))

        #Fpos.extend(Fneg)
        Fpos.extend(Fneg)
        #print "Fpos: " + str(len(Fpos))

        #print "Fpos: " + str(Fpos[:5])
        #print "Fpos: " + str(Fpos[1168:1172])
        #print "Fpos: " + str(Fpos[-5:])


        return Fpos

    def simulateLightcurve(self, fourier, expon=True, lcmean = None):

        ### length of time interval
        tmax = 1.0/self.dnu
        ### sampling rate
        dt = tmax/self.n

        #print(self.n)

        ### make a time array
        time = np.arange(len(fourier))*tmax/self.n


        f = fourier

        phi = np.fft.ifft(f)#/np.sqrt(self.n)#/(self.nphot**0.5)

        phi = np.array([x.real for x in phi])

        ### if expon == True, transform into lognormally distributed 
        ###light curve such that there are no values < 0:
        if expon == False:
            flux = phi 
        elif expon == True and not lcmean == None:
            lncounts = np.exp(phi)
            flux = lncounts*lcmean/np.mean(lncounts)
        else:
            raise Exception("You must either specify a mean flux or set expon=False !")

        lc = lightcurve.Lightcurve(time, counts=flux)

        return lc


    def makePeriodogram(self, fourier, norm='variance'):

       f = fourier 
       f2 = np.array(f).conjugate()
       ff2 = np.real(f*f2)
       s = ff2[0:self.n/2]#*2.0/(self.fnyquist*2*self.fnyquist)

       if norm.lower() in ['variance', 'var']:
           s = s*2.0/(self.fnyquist*2*self.fnyquist)

       if norm.lower() in ['leahy']:
           s = 2.0*s/self.nphot
      
       if norm.lower() in ['rms']:
           s = 2.0*s/(df*self.nphot**2.0)


       ps = powerspectrum.PowerSpectrum()
       ps.freq = self.freq
       ps.ps = s

       return ps

##########################################################
    

    
    






