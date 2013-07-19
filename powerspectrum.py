import numpy as np
import math
import scipy
import scipy.optimize
import lightcurve
import generaltools

# this is the right version of powerspectrum.py

### add up a number of periodograms or average them
## psall: list of periograms (note: *must* have same resolution!)
## method: 
##     'average', 'avg', 'mean': average periodograms
##     'add', 'sum': add them up
def add_ps(psall, method='avg'):

    pssum = np.zeros(len(psall[0].ps))
    for x in psall:
        pssum = pssum + x.ps

    if method.lower() in ['average', 'avg', 'mean']:
        pssum = pssum/len(psall)

    psnew = PowerSpectrum()
    psnew.freq = psall[0].freq
    psnew.ps = pssum
    psnew.n = psall[0].n
    psnew.df = psall[0].df
    psnew.norm = psall[0].norm
    return psnew 

### REDOING THIS IN CLASSES ####
class PowerSpectrum(lightcurve.Lightcurve):
    def __init__(self, lc = None, counts = None, nphot=None, norm='leahy', verbose=False):

        self.norm = norm

        if isinstance(lc, lightcurve.Lightcurve) and counts == None:             
            pass
            # print "You put in an object of type Lightcurve. Go you!"
        elif not lc == None and not counts == None:
            if verbose == True:
                print "You put in a standard light curve (I hope). Converting to object of type Lightcurve"
            lc = lightcurve.Lightcurve(lc, counts, verbose=verbose)
        else:
            self.freq = None
            self.ps = None
            self.df = None
            return
        #    raise Exception("You screwed up. Put in either a Lightcurve or two lists of time, counts")

        if nphot == None:
            nphots = np.sum(lc.counts)
        else:
            nphots = nphot
        nel = np.round(lc.tseg/lc.res)
        print(nel)
        #if lc.tseg/lc.res > nel:
        #    nel = nel + 1
     
        #print "lc.tseg: " + str(lc.tseg)
        #print "lc.res: " + str(lc.res)
        df = 1.0/lc.tseg
        fnyquist = 0.5/(lc.time[1]-lc.time[0])
        #freq = np.arange(nel/2)*df ### make a list of frequencies 
        #print "len(freq): " + str(len(freq))
        #print "len(lc.counts): " + str(len(lc.counts))
        fourier= scipy.fft(lc.counts) ### do Fourier transform
        f2 = fourier.conjugate() ### do conjugate
        ff = f2*fourier   ### multiply both together
        fr = np.array([x.real for x in ff]) ### get out the real part of ff
        ps = 2.0*fr[0: (nel/2 )]/nphots
        freq = np.arange(len(ps))*df
        #print "ps0 : " + str(ps[0])
        #print "df : " + str(df)
        #print "nel: " + str(nel)
        #print "nphots" + str(nphots)
 

        if norm.lower() in ['leahy']:
            self.ps = ps#2.0*fr[0: (nel/2 -1)]/nphots   ### do the power spectrum
            
        elif norm.lower() in ['rms']:
            #self.ps = 2.0*lc.tseg*2.0*np.sqrt(nel)*fr[0:(nel/2 -1)]/(nphots**2.0)
            self.ps = ps/(df*nphots)        
            print "ps rms: " + str(self.ps[0])
        ### put frequency to mid-frequencies
        elif norm.lower() in ['variance', 'var']:
            #self.ps = ps*nphots/len(lc.counts)
            self.ps = ps*nphots/len(lc.counts)**2.0

            #self.ps = ps*nphots/(2.0*fnyquist)**2.0


        self.freq = [f+(freq[1]-freq[0])/2.0 for f in freq]
        self.df = self.freq[1] - self.freq[0]
        self.nphots = nphots
        self.n = len(lc.counts)

    def rebinps(self, res, verbose=False):
        ### frequency range of power spectrum
        flen = (self.freq[-1] - self.freq[0])
        ### calculate number of new bins in rebinned spectrum
        bins = np.floor(flen/res) 
        ### calculate *actual* new resolution
        self.bindf = flen/bins
        ### rebin power spectrum to new resolution
        binfreq, binps, dt = self._rebin(self.freq, self.ps, bins, method='mean', verbose=verbose)
        newps = PowerSpectrum()
        newps.freq = binfreq
        newps.ps = binps
        newps.df = dt

        return newps


    ### NEED TO TEST THIS!! DON'T KNOW WHETHER IT WORKS!!!
    def rebingeo(self, factor, method = 'sum'):
        freq = np.array(self.freq)
        ps = np.array(self.ps)

        ### initial frequency resolution
        df = freq[1]-freq[0]

        ###low index boundaries
        ilow = 1
        ihigh = 2
        #### new frequency array, keep first frequency (N_phot) unchanged
        freqgeo = [freq[0]]
        psgeo = [ps[0]]

        ### start at old bin boundary of first bin
        dffrac = 0.0
        fracbin = 0.0

        dfnew = df
        dfnewlist = []
        #print "dfnew before loop: " + str(dfnew)
        #print "Factor: " + str(factor)
        while ihigh < len(freq)-1 and ilow < len(freq)-1:
            ### new frequency spacing for next bin
            #dfnew = dfnew*(1.0+factor)
        #    print "dfnew: " + str(dfnew)
        #    print "ilow: " + str(ilow)
            ### part to add to index if fracind is not zero
            dftemp = dfnew - dffrac
            ### fraction of old lower bin to be added to new bin
            fracbinlow = 1.0 - fracbin
            freqbin = freq[ilow+1] + 0.5*dftemp
            ### index after which the frequency with freq[ilow] + dftemp would be sorted in
            freqmax = freq[ilow] + dftemp
        #    print "freqmax: " + str(freqmax)
            ihigh = freq.searchsorted(freq[ilow]+dftemp)
            if ihigh >= len(freq):
                ihigh = len(freq) - 1
        #    print "ihigh: " + str(ihigh)
            ### fraction of last bin
            dffrac = freqbin - freq[ihigh]
            ### fraction of old bin that is in new bin
            fracbin = dffrac/df
            dfnew = dfnew*(1.0+factor)

            ### add together powers
            ### is this right, or should it be ihigh+1 ????
            print "method: " + str(method)
            if method == 'sum':
                psnew = sum(ps[ilow+1:ihigh]) + ps[ilow]*fracbinlow + ps[ihigh]*fracbin 
            elif method.lower() in ['average', 'mean']:
                pscomp = ps[ilow+1:ihigh]
                list(pscomp).append(ps[ilow]*fracbinlow)
                list(pscomp).append(ps[ihigh]*fracbin)
                psnew = np.average(pscomp)
            freqgeo.append(freqbin)
            psgeo.append(psnew)
            dfnewlist.append(dfnew) 
            ilow = ihigh

                
        return freqgeo, psgeo, dfnewlist
    ##########################################


    def findmaxpower(self):
        psfiltered = filter(lambda x: x >= 100.0, self.ps)
        maxpow = max(psfiltered)
        return maxpow

    def fitps(self):
        freq = self.freq
        ps = self.ps

        noiseminfreq = np.array(freq).searchsorted(1000.0)
        noise = np.mean(self.ps[noiseminfreq:])

        print "The noise level above 1000 Hz is: " + str(noise)
        if noise < 1.95 or noise > 2.05:
            self.deadtime = True
        else:
            self.deadtime = False

        ### fit just low-frequency power law
        func = lambda x,a,b: a*x+b

        ### pick out only frequencies < 100 Hz for the power law 
        ### at least until I find a way to characterize the break in the spectrum
        maxpl = np.array(freq).searchsorted(100.0)
        
        logf = np.array([math.log10(f) for f in freq])
        logps = np.array([math.log10(p) for p in ps])

        ### guesses for initial values
        p_init = [-2.0, 10000.0]

        ### do non-linear least-squares fitting
        ### covar contains the estimated covariance of pfinal
        ### The diagonals provide the variance of the parameter estimate
        if verbose == True:
            print "Fitting power law  + constant with initial guesses " + str(p_init) + " ..."
        pfinal, covar = scipy.optimize.curve_fit(func, logf[1:maxpl+1], logps[1:maxpl+1], p0=p_init)
        if verbose == True:
            print("The best-fit parameters for a function of type b*x**a + c are: ")
            print("index = " + str(pfinal[0]))
            print("normalization = " + str(pfinal[1]))
 
        ### fitted power law index and normalization
        pind = pfinal[0]
        pnorm = 10**pfinal[1]

        ### make power law along whole spectrum
        ypl = 10.0**pfinal[1]*(np.array(freq)**(pfinal[0]))

        ### reverse power law to be able to use searchsorted in the next step
        ypl = list(ypl)
        ypl.reverse()
        
        ### find the frequency where powerlaw and noise level meet
        plstop = np.array(ypl).searchsorted(noise)

        ### reverse power law back into proper order
        ypl.reverse()

        maxpl = len(freq) - plstop
        print("Frequency where power law and noise are equal: " + str(freq[maxpl]) + " Hz.")
        yfit = ypl[:maxpl]
        yfit.extend([noise for x in freq[maxpl:]])
        self.psfit = yfit
        self.powerlaw = ypl
        self.noiselevel=noise

        return pind, pnorm



    def variability(self, verbose=False):
        ### function to be fit:
        ### combination of power law and constant (should be 2, but not neccessarily)
        #func = lambda x, a, b, c: b*(x**a) + c
 
        freq = self.freq
        ps = self.ps

        noiseminfreq = np.array(freq).searchsorted(1000.0)
        noise = np.mean(self.ps[noiseminfreq:])
        print "The noise level above 1000 Hz is: " + str(noise)
        #func = lambda x, a, b: b*(x**a) + noise
        func = lambda x,a,b: a*x+b
        
        ### find maximum power for power law:
        maxpl = np.array(freq).searchsorted(100.0)
        logf = np.array([math.log10(f) for f in freq])
        logps = np.array([math.log10(p) for p in ps])      
 
        ### guesses for initial values
        p_init = [-2.0, 10000.0]

        ### do non-linear least-squares fitting
        ### covar contains the estimated covariance of pfinal
        ### The diagonals provide the variance of the parameter estimate
        if verbose == True:
            print "Fitting power law  + constant with initial guesses " + str(p_init) + " ..."
        pfinal, covar = scipy.optimize.curve_fit(func, logf[1:maxpl+1], logps[1:maxpl+1], p0=p_init)
        if verbose == True:
            print "The best-fit parameters for a function of type b*x**a + c are: "
            print "a = " + str(pfinal[0])
            print "b = " + str(pfinal[1])
            #print "c = " + str(pfinal[2])

            #print "Note that c is the average noise level."
            #print "Deviation from Poisson noise by " + str(2.0-pfinal[2])

        ### just the fitted power law
        pind = pfinal[0]
        pnorm = 10**pfinal[1]

        ypl = 10.0**pfinal[1]*(np.array(freq)**(pfinal[0]))
        print ypl
        #ypl = 10.0**ypl
        #print ypl
        ypl = list(ypl)
        ypl.reverse()
        print "ypl: " + str(ypl)
        plstop = np.array(ypl).searchsorted(noise)
        ypl.reverse()
        print "plstop: " + str(plstop)
        maxpl = len(freq) - plstop
        print "maxpl: " + str(maxpl)
        yfit = ypl[:maxpl]
        yfit.extend([noise for x in freq[maxpl:]])
        #yconst = [pfinal[2] for bla in freq]
        ### complete fit
        #yfit = ypl + pfinal[2]


        ### from the power law, find the frequency where the power low drops below the
        ### constant noise level
        ### list of all indices of values in ypl that are smaller than the const. noise level
        #lowvalues = [list(ypl).index(y) for y in ypl if y < pfinal[2]]

        #params = {'backend':'ps', 'text.fontsize':10, 'legend.fontsize':10}
        #plt.rcParams.update(params)

        #print "Done plotting ..."


        ### *first* value where power law dips below noise level
        #varts_index = lowvalues[0]
        ### corresponding frequency:
        #varfreq = freq[varts_index]

        #print "The variability time scale is F = " + str(varfreq) + " Hz, t = " + str(1.0/varfreq) + " s."

        self.psfit = yfit
        self.powerlaw = ypl
        self.noiselevel=noise

        return pind, pnorm
        #return varfreq

    def checknormal(self, freq, ps):
        ### checks the normalization of a power spectrum above fnyquist/10 Hz
        fmin = max(freq)/10.0
        minind = np.array(freq).searchsorted(fmin)
        psnew = ps[minind:-1]
        normlevel = np.average(psnew)
        normvar = np.var(psnew)

        return normlevel, normvar




#### AUXILIARY FUNCTIONS #######

#def errparam()



