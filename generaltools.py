############################################################################
# GENERAL TOOLS USEFUL OR NECESSARY FOR DATA PROCESSING
# AND ANAYLSIS
#
# last changed: 
#  - July 24, 2013: added function choice_hack as a substitute for np.choice
#  - March 24, 2012: Added documentation
#  - October 24, 2011: worked out kinks in barycentering
# 
##################

#### IMPORTS #################
#
## STANDARD LIBRARIES
# These are standard libraries built
# into python 2.7
# --> if you don't have python 2.7, some of these might not work!

from __future__ import with_statement
from collections import defaultdict
import argparse ### new in python 2.6 (I think)
import cPickle as pickle
import fractions
import math
import os
import glob


#### COMMON PACKAGES 
#
# Commonly installed libraries like
# numpy, scipy or pyfits
#
# Installatiion required except for matplotlib
#

### find out if matplotlib is installed. 
### If not, don't do any plots.
try:
#    import matplotlib
#    matplotlib.use('PS')
    import matplotlib.pyplot as plt
    from pylab import *
except ImportError: 
    matplotlib = None

import numpy as np
import scipy


#### OWN MODULES
# 
# Modules I wrote myself,
# like this one :-) 
#
#

### JPL Ephemeris scripts
try:
    import jplephread
    import jplephinterp
except ImportError:
    print("You don't have jplephread and jplephinterp! You can't do barycentering like this!")

### AUTOCORRELATION FUNCTION ####
#
#  Computes the autocorrelation function, i.e. the correlation
#  of a data set with itself.
#  To do this, shift data set by one bin each time and compute correlation for
#  the data set with itself, shifted by i bins
#
#  If the data is _not_ correlated, then the autocorrelation function is the delta
#  function at lag = 0 
#
#  The autocorrelation function can be computed explicitly, or it can be computed
#  via the Fourier transform (via the Wiener-Kinchin theorem, I think) 
#
#  VARIABLES:  
#
#  x [list]: data set
#  nlags [int]: maximum lag, i.e. number of times by which the data set should maximally be shifted
#  fourier [bool]: compute ACF via Fourier transform (True) or explicitly (False)
#  norm [bool]: normalize ACF to 1
#
def autocorr(x,nlags = 100, fourier=False, norm = True):

    ### empty list for the ACF
    r = []
    ### length of the data set
    xlen = len(x)

    ### shift data set to a mean=0 (otherwise it comes out wrong) 
    x1 = np.copy(x) - np.mean(x)
    x1 = list(x1)

    ### add xlen zeros to the array of the second time series (to be able to shift it)
    x1.extend(np.zeros(xlen))

    ### if not fourier == True, compute explicitly
    if not fourier:
        ### loop over all lags
        for a in range(nlags):
            ### make a np.array of 2*xlen zeros to store the data set in
            x2 = np.zeros(len(x1))
            ### put data set in list, starting at lag a
            x2[a:a+xlen] = x-np.mean(x)
            ### compute autocorrelation function for a, append to list r
            r.append(sum(x1*x2)/((xlen - a)*np.var(x)))      
         
    ### else compute autocorrelation via Fourier transform
    else:
        ### Fourier transform of time series
        fourier = scipy.fft(x-np.mean(x))
        ### take conjugate of Fourier transform
        f2 = fourier.conjugate()
        ### multiply both together to get the power spectral density
        ff = f2*fourier
        ### extract real part
        fr = np.array([b.real for b in ff])
        ps = fr
        ### autocorrelation function is the inverse Fourier transform
        ### of the power spectral density
        r = scipy.ifft(ps)
        r = r[:nlags+1]
    ### if norm == True, normalize everything to 1
    if norm:
        rnew = r/(max(r))
    else:
        rnew = r
    return rnew  

#### Correlation function
#
# NOT SURE THIS IS CORRECT! CHECK THIS 
# WHEN I HAVE ACCESS TO INTERNET!
#
def correlate(x, y, nlags=100, norm=True):

    ### empty list for the ACF
    r = []
    ### length of the data set
    xlen = len(x)

    ### shift data set to a mean=0 (otherwise it comes out wrong) 
#    x1 = np.copy(x) - np.mean(x)
    x1 = y - np.mean(y)
    x1 = list(x1)

    ### add xlen zeros to the array of the second time series (to be able to shift it)
    x1.extend(np.zeros(xlen))

    ### if not fourier == True, compute explicitly
    for a in range(nlags):
            ### make a np.array of 2*xlen zeros to store the data set in
            x2 = np.zeros(len(x1))
            ### put data set in list, starting at lag a
            x2[a:a+xlen] = x-np.mean(x)
            ### compute autocorrelation function for a, append to list r
            r.append(sum(x1*x2)/((xlen - a)*np.var(x)))
    ### if norm == True, normalize everything to 1
    if norm:
        rnew = r/(max(r))
    else:
        rnew = r
    return rnew


###### QUICK AND DIRTY REBINNING OF LIGHT CURVES #####################
#
#
# Does a quick and dirty rebin of light curves by integer numbers.
#
#
#
#
#
#
#
def rebin_lightcurve(times, counts, n=10, type='average'):

    nbins = int(len(times)/n)
    dt = times[1] - times[0]
    T = times[-1] - times[0] + dt
    bin_dt = dt*n
    bintimes = np.arange(nbins)*bin_dt + bin_dt/2.0 + times[0]

    nbins_new = int(len(counts)/n)
    counts_new = counts[:nbins_new*n]
    bincounts = np.reshape(np.array(counts_new), (nbins_new, n))
    bincounts = np.sum(bincounts, axis=1)
    if type in ["average", "mean"]:
        bincounts = bincounts/np.float(n)
    else:
        bincounts = bincounts

    #bincounts = np.array([np.sum(counts[i*n:i*n+n]) for i in range(nbins)])/np.float(n)
    #print("len(bintimes): " + str(len(bintimes)))
    #print("len(bincounts: " + str(len(bincounts)))
    if len(bintimes) < len(bincounts):
        bincounts = bincounts[:len(bintimes)]

    return bintimes, bincounts



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

##############################################################



#### APPROXIMATE HESSIAN ######################################3
#
# Note: I stole this routine from 
# 	scikits.statsmodels.sandbox.regression.numdiff
#       and should give credit where credit is due ...
#
#
#
def approx_hess3(x0, f, epsilon=None, args=(), **kwargs):
    '''calculate Hessian with finite difference derivative approximation

    Parameters
    ----------
    x0 : array_like
       value at which function derivative is evaluated
    f : function
       function of one array f(x)
    epsilon : float
       stepsize, if None, then stepsize is automatically chosen

    Returns
    -------
    hess : ndarray
       array of partial second derivatives, Hessian

    Notes
    -----
    based on equation 9 in
    M. S. RIDOUT: Statistical Applications of the Complex-step Method
    of Numerical Differentiation, University of Kent, Canterbury, Kent, U.K.

    The stepsize is the same for the complex and the finite difference part.
    '''
    EPS = 2.2204460492503131e-016


    if epsilon is None:
        h = EPS**(1/5.)*np.maximum(np.abs(x0),1e-2) # 1/4 from ...
    else:
        h = epsilon
    xh = x0 + h
    h = xh - x0
    ee = np.diag(h)
    hess = np.outer(h,h)

    n = dim = np.size(x0) #TODO: What's the assumption on the shape here?

    if "neg" in kwargs.keys():
        for i in range(n):
            for j in range(i,n):
                hess[i,j] = (f(*((x0 + ee[i,:] + ee[j,:],)+args), neg=kwargs["neg"][0])
                                - f(*((x0 + ee[i,:] - ee[j,:],)+args), neg=kwargs["neg"][0])
                             - (f(*((x0 - ee[i,:] + ee[j,:],)+args), neg=kwargs["neg"][0] )
                                - f(*((x0 - ee[i,:] - ee[j,:],)+args), neg=kwargs["neg"][0]))
                             )/4./hess[i,j]
                hess[j,i] = hess[i,j]
    else:
        for i in range(n):
            for j in range(i,n):
                hess[i,j] = (f(*((x0 + ee[i,:] + ee[j,:],)+args))
                                - f(*((x0 + ee[i,:] - ee[j,:],)+args))
                             - (f(*((x0 - ee[i,:] + ee[j,:],)+args))
                                - f(*((x0 - ee[i,:] - ee[j,:],)+args)))
                             )/4./hess[i,j]
                hess[j,i] = hess[i,j]



    return hess






#### PHOTON CLASS ###############################################
#
# General Photon class, with subclasses for specific instruments.
# Instances of these class are objects with at least two attributes:
# - photon arrival time
# - energy/ channel the photon was detected in
#
# Subclasses: gbm.GBMPhoton, rxte.RXTEPhoton, SatPos
#
#
#
class Photon(object):
    def __init__(self, time, energy):
         self.time = time
         self.energy=energy

    ### convert mission time to Modified Julian Date
    def mission2mjd(self, mjdrefi, mjdreff, timezero=0.0):
        self.mjd = (mjdrefi + mjdreff) + (self.time + timezero)/86400.0

    ### Auxiliary function for computation of energy boundaries
    ### DON'T CALL ON OBJECT, CALL WITHIN APPROPRIATE METHOD!
    def _in_range(self, lower, upper):
        if lower <= self.energy <= upper:
            return True
        else:
            return False

###################################################################


def _checkinput(gti):
    if len(gti) == 2:
        try:
            iter(gti[0])
        except TypeError:
            return [gti]
    return gti



#### MAKE A GENERAL DATA OBJECT ###################################
#
# This is a general class for X-ray and Gamma-Ray data with methods
# for filtering data.
# 
# Note: Strictly speaking, this serves as a superclass for its subclasses,
# supplying common attributes and methods.
# DON'T CALL THIS CLASS BY ITSELF! RATHER CALL A SUBCLASS!
#
# Subclasses:	 gbm.GBMData, rxte.RXTEData 
#
#
#
class Data(object):
    def __init__(self):
        raise Exception("Don't run this! Use subclass RXTEData or GBMData instead!")

    ### Filter out photons that are outside energy thresholds cmin and cmax
    def filterenergy(self, cmin, cmax):
        self.photons= [s for s in self.photons if s._in_range(cmin, cmax)]

    ### For instruments supplying good time intervals, this function enables 
    ### GTI filtering
    ### GTIs can eitehr be passed into the function or are an attribute of the
    ### data subclass used.
    def filtergti(self, gti=None):
        if not gti:
            gti = self.gti
        gti=_checkinput(gti)
        #print gti
        filteredphotons = []
        ### Use _unbarycentered_ time to filter GTIs
        ### NEED TO CHECK WHETHER THAT IS TRUE FOR BOTH
        ### RXTE AND GBM !!! 
        times = np.array([t.unbarycentered for t in self.photons])
        ### note: this method below is more clunky than using filter(),
        ### but it's much faster, too! :-)
        for g in gti:
            tmin = times.searchsorted(g[0])
            tmax = times.searchsorted(g[1])
            photons = self.photons[tmin:tmax]
            filteredphotons.extend(photons)
        self.photons = filteredphotons


    ### NEED TO TEST THIS PROPERLY!
    ### Filter for a burst, using a tuple or list in bursttimes
    ### if blen is not given, then it is calculated from bursttimes.
    ### the flag 'bary' sets whether bursttimes are barycentered.
    def filterburst(self, bursttimes, blen=None, bary=False):
        tstart= bursttimes[0]
        tend = bursttimes[1]
        if blen == None:
            blen = tend - tstart

        #tunbary = np.array([s.unbarycentered for s in self.photons])
        time = np.array([s.time for s in self.photons])

        ### The bary flag sets whether the burst times are barycentered
        ### or not. By default (especially for RXTE data), this is not the
        ### case. For GBM data, it usually is
        if bary == False:
            tunbary = np.array([s.time for s in self.photons])
            stind = tunbary.searchsorted(tstart)
            eind = tunbary.searchsorted(tend)
        else:
            stind = time.searchsorted(tstart)
            eind = time.searchsorted(tend)

        self.burstphot = self.photons[stind:eind]

    ### Barycenter Times of Arrival ###
    ### Note that this needs barycentered position history data
    ### supplied in the form of a PosHist object (or relevant subclasses)
    def obsbary(self, poshist):

        ### photon times of arrival in MET seconds
        tobs = np.array([s.time for s in self.photons])
        ### barycentered satellite position history time stamps in MET seconds
        phtime = np.array([s.phtime for s in poshist.satpos])
        ### Interpolate barycentering correction to photon TOAs
        tcorr = np.interp(tobs, phtime, poshist.tdiff)
        ### barycentered photon TOAs in TDB seconds since MJDREFI
        ctbarytime = (tobs + tcorr)
        ctbarymet = ctbarytime + self.mjdreff*8.64e4  # barycentered TOA in MET seconds

        ### barycenter trigger time the same way as TOAs
        trigcorr = np.interp(self.trigtime, phtime, tdiff)
        trigcorr = (self.trigtime + trigcorr)

        ### barycentered photon TOAs in seconds since trigger time
        ctbarytrig = ctbarytime - trigcorr 
        ### barycentered photon TOAs as Julian Dates
        ctbaryjd = ctbarytime/8.64e4 + self.mjdrefi + 2400000.5

        ### return dictionary with TOAs in different formats
        ctbary = {'barys': ctbarytime, 'trigcorr': trigcorr, 'barymet': ctbarymet, 'barytrig': ctbarytrig, 'baryjd': ctbaryjd}
        return ctbary


######################################################################
######################################################################

##### SATELLITE POSITION CLASS ###############################
#
# This is a very simple subclass of Photon that defines a Satellite Position
# Object with four attributes: a time stamp and three spatial
# dimensions.
# To be used within the framework of the PosHist class.
# Can be extended should there be demand.
#
class SatPos(Photon):
    def __init__(self, phtime,time, x, y, z):
        self.time = time
        self.phtime = phtime
        self.x = x
        self.y = y
        self.z = z
###############################################################




##### SATELLITE POSITION *HISTORY* CLASS #####################
#
# Superclasses: -Data
#
# Objects of this class get data from a Position History File 
# (depending on satellite) and stores it as list of SatPos objects.
# Barycentering capability included.
# Barycentering needs Ra and Dec as well as fbldata and 
# JPLEPH.200 on disk, preferably in the same folder
# as the scripts.
#
#
class PosHist(Data):

    ### initialize Position History object:
    ### need: - A Position History file name
    ###       - Ra and Dec of object (for barycentering)
    ###       - the fbldata in a dictionary *OR*
    ###       - path for fbldata (without filename) if it's not in the same directory
    ### Note: You probably don't want to run this, but use 
    ### a subclass instead.

    def __init__(self, ra = None, dec = None, fbldata = None, fblpath = None):

        ### store RA and Dec for barycentering
        if not ra == None and not dec == None:
            self.ra = ra
            self.dec = dec
            self.bary = True
        else:
            print "No RA and Dec given. Will not barycenter!"
            self.bary = False

        if not fbldata == None:
            PosHist.fbldata = fbldata
        else:
            PosHist.fbldata = self._readfbl(fblpath=fblpath)

        if self.bary == True:
            try:
                print "Barycentering satellite position history ..."
                self.tdiff = self.baryposhist(PosHist.fbldata)
            except AttributeError:
                raise Exception("You don't want to run this! Run subclass instead!")
                      

    ### Auxiliary method to read in fbldata from file.
    ### Note: Don't run this method on object, but use in 
    ### dedicated method that needs it! 
    def _readfbl(self, fblpath = None):

        if fblpath == None:
            fbldir = os.path.dirname(__file__)
        else:
            fbldir = fblpath
        fblname = os.path.join(fbldir, "fbldata.dat")
        f=open(fblname, 'r')
        const0, freq0, phase0, iterms= [], [], [], []
        for line in f:
            if line.startswith('#'):
                iterms.append(len(const0))
            else:
                number1, number2, number3 = [float(x[:-1]) for x in line.split()]
                const0.append(number1)
                freq0.append(number2)
                phase0.append(number3)
        texp=np.zeros(len(const0))
        texp[iterms[0]:iterms[1]]=1.0
        texp[iterms[1]:iterms[2]]=2.0
        texp[iterms[2]:iterms[3]]=3.0
        texp[iterms[3]:]=4.0
        output_lists={'const0':np.array([const0]), 'freq0':np.array([freq0]), 'phase0':np.array([phase0]), 'iterms':iterms, 'texp':texp}
        return output_lists


    ### THIS BARYCENTERS POSITION HISTORY DATA, BUT IT COULD IN THEORY BARYCENTER OTHER
    ### POSITION HISTORY DATA THAN FERMI GBM'S!
    def baryposhist(self, fbldata):

        ### read ephemeris coefficients from file
        ### as default, JPLEPH.200 is used and it is assumed that file
        ### resides in the same directory as the python script.
        ### It is possible to use another file/directory if specified,
        ### but I have no experience with any other data than JPLEPH.200, so you're on
        ### your own there :-) 
        info, coeffs = jplephread.jplephread(self.jdlimits, ephname = 'JPLEPH.200', ephdir = 'current')

        ### speed of light in km/s (I think)
        c=info['c']/1.0e3

        ### position history time stamps, adding MJDREFF in seconds
        ### to get time in seconds since MJDREFI (=base_jd)
        time = np.array([s.phtime for s in self.satpos]) + self.sec_offset

        ### position of space craft at time stamps in time
        position = np.array([[s.x, s.y, s.z] for s in self.satpos])

        ### Base JD = 2400000.5 + MJDREFI
        tbase = self.basejd

        ### test whether RA and Dec are given, complain if not!
        try:
            rarad = math.pi*float(self.ra)/180.0
            decrad = math.pi*float(self.dec)/180.0
        except AttributeError:
            raise Exception('You forgot to give Ra and Dec! Try again!')


        ### list with direction of the source in three dimensions
        srcdir = [math.cos(decrad)*math.cos(rarad) , math.cos(decrad)*math.sin(rarad), math.sin(decrad)]

        ### tdb 2 tt conversion
        dtein = self._tdb2tdt(time/8.64e4, tbase, fbldata)

        scnew = [(a+b)/8.64e4 for a,b in zip(time, dtein)] ### this is in JD since base_jd

        ### earth position and velocity
        xe, ye, ze, vxe, vye, vze = jplephinterp.jplinterp(3, info, scnew, tbase, coeffs)

        ### transform to barycentric reference frame from geocentric frame
        dtlorentz = (vxe*position[:,0] + vye*position[:,1] + vze*position[:,2])/(1.0e3*c**2.0)

        ### transformed times in seconds since MJDREFI
        sctdb = time + dtlorentz + dtein 

        ### transformed time stamps in JD since MJDREFI
        sctdbjd = sctdb/8.64e4

        ### position (and velocities, although they aren't needed) of the sun
        xs, ys, zs, vxs, vys, vzs = jplephinterp.jplinterp(11, info, sctdbjd, tbase, coeffs)

        xsc, ysc, zsc, delx, dely, delz = [], [], [], [],[],[]
        for i in range(len(xe)):
            xsc.append(xe[i] + position[i,0]/1.0e3)
            ysc.append(ye[i] + position[i,1]/1.0e3)
            zsc.append(ze[i] + position[i,2]/1.0e3)
            delx.append(xsc[i] - xs[i])
            dely.append(ysc[i] - xs[i])
            delz.append(zsc[i] - xs[i])

        xsc = np.array(xsc)
        ysc = np.array(ysc)
        zsc = np.array(zsc)
     
        delx = np.array(delx)
        dely = np.array(dely)
        delz = np.array(delz)
  
        ### geometric correction 
        dtgeo = ((srcdir[0]*xsc + srcdir[1]*ysc + srcdir[2]*zsc)/c)

        #delx = xsc - xs
        #dely = ysc - xs
        #delz = zsc - xs

        ### calculate Shapiro delay
        costheta = (srcdir[0]*delx + srcdir[1]*dely + srcdir[2]*delz)/(delx**2 + dely**2 + delz**2)**0.5
        dtshapiro = -9.8509819e-6 * np.log(1+costheta)

        ### barycentered position history times
        self.barytdb = sctdb+ dtgeo - dtshapiro
        ### and difference between barycentered and unbarycentered times
        tdiff = self.barytdb - np.array([s.phtime for s in self.satpos])
        return tdiff


    ### Auxiliary method converting between tdb and tt
    def _tdb2tdt(self, time, tbase, fbldata):
        t = ((tbase - 2451545.0) + time)/365250.0
        const = fbldata['const0']
        freq= fbldata['freq0']
        phase = fbldata['phase0']
        texp = fbldata['texp']
        dt = []
        for i, tn in enumerate(t):
            ph = freq*tn + phase
            sint = np.sin(ph)
            bla = const*(tn**(texp))
            quot = bla*sint
            dt.append(np.sum(quot)*1.0e-6)
        return dt


######################################################################




####### DEFINE SOME UNIVERSAL COMMAND LINE OPTIONS ##################
#
#  !!! UNDER CONSTRUCTION !!!
#

arglibrary = {"-f": dict(action="store", dest="filename", help='Enter a file name to be used.'), "-o": dict(action="store", dest="outfile", help='Enter the filename for the output file.'), "-i": dict(action="store", dest="obsid", help='Enter the burst number (obsid) to be used'), "-d": dict(action="store", dest="date", help='Enter a date for the position history file. Format: DD/MM/YYYY'), "--det": dict(action="store", dest="detfile", help='Enter a filename for the file with detector choices.')}
arglibrary["--filename"] = arglibrary["-f"]
arglibrary["--outfile"] = arglibrary["-o"]
arglibrary["--obsid"] = arglibrary["-i"]
arglibrary["--date"] = arglibrary["-d"]


def add_arguments(parser, *arglist):
    for arg in arglist:
        parser.add_argument(arg, **arglibrary[arg])
    return parser

#parent_parser = argparse.ArgumentParser(add_help=False)
#parent_parser.add_argument('-f', '--filename', action="store", dest="filename", help='Enter a filename to be used.')
#parent_parser.add_argument('-o', '--outfile', action="store", dest="outfile", help='Enter the filename of the resulting output file.')
#parent_parser.add_argument('-i', '--obsid', action="store", dest="obsid", help = "Enter the burst number (obsid) here.")
############################################################################



#### READ ASCII DATA FROM FILE #############
#
# This is a useful little function that reads
# data from file and stores it in a dictionary
# with as many lists as the file had columns.
# The dictionary has the following architecture
# (e.g. for a file with three columns):
#
# {'0':[1st column data], '1':[2nd column data], '2':[3rd column data]}
#
#
# NOTE: Each element of the lists is still a *STRING*, because
# the function doesn't make an assumption about what type of data
# you're trying to read! Numbers need to be converted before using them!
#
def conversion(filename):
    f=open(filename, 'r')
    output_lists=defaultdict(list)
    for line in f:
        if not line.startswith('#'):
             line=[value for value in line.split()]
             for col, data in enumerate(line):
                 output_lists[col].append(data)
    return output_lists




# ASK FOR NUMERIC VALUE ################
def ask_input(prompt):
    input=raw_input(prompt)
    return input
########################################

# THIS ASKS WHETHER THE USER WANTS TO SPECIFY A VALUE #####
def ask_limits(prompt, retries=3, complaint='yes or no, please!'):
    while True:
          ok=raw_input(prompt)
          if ok in ('y', 'ye', 'yes'):
              input=ask_input('Please specify an input value: ')
              return input
          if ok in ('n', 'no', 'nop', 'nope'):
              return False
          if retries<0 : raise IOError, 'refusenik user'
          print complaint
######################################################################


#### MAKE A SIMPLE PLOT
#
# REQUIRES:     - matplotlib
#               - pyplot
#               - postscript
#
def makeplot(x, y, xlabel, ylabel, title, xrange, yrange, log, outfile):
    if matplotlib == None:
        print "Matplotlib not installed. Returning ..."
        return
    plt.plot(x, y, color='black')
    if log == True:
        plt.yscale('log')
        plt.xscale('log')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.axis([xrange[0], xrange[1], yrange[0], yrange[1]])
    plt.title(title)
    plt.savefig(outfile)
    plt.close()
    return

##### GET DATA FROM PICKLED PYTHON OBJECT (FROM PROCESSING PIPELINE)
#
# Pickling data is a really easy way to store data in binary format.
# This function reads back pickled data and stores it in memory.
#
#
def getpickle(picklefile):
    file = open(picklefile, 'r')
    procdata = pickle.load(file)
    return procdata
########################################################################



#######################################################################
#
#
#
#
#
#
#
class TwoPrint(object):

    def __init__(self,filename):
        self.file = open(filename, "w")
        self.filename = filename
        self.file.write("##\n")
        self.close()
        return

    def __call__(self, printstr):
        print(printstr)
        self.file = open(self.filename, "a")
        self.file.write(printstr + "\n")
        self.close()
        return

    def close(self):
        self.file.close()
        return









##### FIND LOWEST COMMON MULTIPLE ####
#
# Find lowest common multiple (LCM) by way of
# findign the greatest common divisor (GCD).
#
# I will need this for a script that will worry about
# artificial signals in smoothed light curves.
#

# Calculate greatest common divisor using Euclid's method
def gcd(a,b):
    while b:
        a,b = b, a % b
    return a

# Now use this to calculate lowest common multiple
def lcm(a,b):
    return a * b // gcd(a, b)


######################################################
#
# NOTE: THIS IS PART OF CLASS LIGHT CURVE
# I ONLY NEED IT UNTIL MY SCRIPTS ARE STILL CHANGED
# TO A CLASS-BASED IMPLEMENTATION
#
# !!! LEGACY CODE! DON'T USE! !!!
#
#def rebin(time, counts, nbins, method = 'sum'):


    ### nbins is the number of bins in the new light curve
#    nbins = int(nbins)
#    #print "nbins: " + str(nbins)
##    print "time: " + str(len(time))
#
    ### told is the _old_ time resolution
#    told = time[1] - time[0]

    ### tseg: length of the entire segment
#    tseg = time[-1] - time[0] + told

    ### move times to _beginning_ of each bin
#    btime = np.array(time) - told/2.0

    ### dt: new time resolution
#    dt = float(tseg)/float(nbins)

    ### check whether old time resolution is larger than new time resolution
#    if dt <= told:
#        print "Old time resolution bigger than new time resolution."
#        print "That's not implemented yet. Returning power spectrum with original resolution."
#        return time, counts, told

    ### tnew is the ratio of new to old bins
#    tnew = dt/told

    ### new array with bin midtimes
#    bintime = [0.5*dt + t*dt for t in range(nbins)]

    ### this fraction is useful, because I can define a number of old bins until the 
    ### boundaries of the old and new bins match up again
    ### this makes it easier in some cases
#    tnewfrac = fractions.Fraction(tnew)

#    top = tnewfrac.numerator
#    bottom = tnewfrac.denominator

    ### if the fraction turns out insanely big (happens due to rounding errors), then I do
    ### only one iteration (since tseg/tseg = 1)
#    if top > tseg:
#        top = tseg
#        bottom = nbins
    #print "top: " + str(top)
    #print "bottom: " + str(bottom)

#    cbin = []

    ### now iterate over all cycles
    #print "tseg/top: " + str(len(tseg/top))



#    for i in range(int(tseg/top)):

        ### I need this index to remember where I left off during the iteration
#        before_ind = 0
        # print "i: " + str(i)
        ### for each cycle iterate through the number of new bins in that cycle
#        for j in range(bottom):
            # print "j: " + str(j)
            ### in the first round, start at the lower edge of the bin:
#            if before_ind == 0:
                #print "tnew: " + str(tnew)
                ### this is the first index to use
#                i0 = int(i*top)
                #print "i0: " + str(i0)
                ### first I sum up all complete old bins in that new bin
#                aint = sum(counts[i0:int(i0+math.floor(tnew))])
                #print "lower index: " + str(i0)
                #print "upper index: " + str(int(i0+math.floor(tnew)))
                #print "values to sum: " + str(counts[i0:int(i0+math.floor(tnew))])

                ### then I set the index of the old bin that is _not_ completely in the new bin
#                fracind = int(i0 + math.floor(tnew) )
                #print "fracind 1 : " + str(fracind)

                ### frac is the fraction of the old bin that's in the new bin
#                frac = tnew - math.floor(tnew)
                #print "tnew fractional part: "  + str(tnew- math.floor(tnew))

                ### if frac is not zero, then compute fraction of counts that goes into my new bin
#                if frac < 1.0e-10:
#                    frac =0

#                if not frac == 0:
#                    afrac = frac*counts[fracind]
                    #print "afrac: " + str(afrac)
#                    cbin.append(aint + afrac) ### append to array with new counts
#                else:
#                    cbin.append(aint)
                #print "cbin: " + str(cbin[-1])

                ### reset before_ind for next iteration in j
#                before_ind = fracind
                #print "before_ind 1 : " + str(before_ind)
#            else:

                ### This new bin doesn't start at the same position as the old bin, hence I start with the fraction
                ### afrac1 is the rest of the preceding old bin that was split up
#                afrac1 = (1.0 - frac)*counts[before_ind]
                # print "afrac1: " + str(afrac1)
                ### 1.0-frac of the bin already done, so define new length for the rest: ttemp 
#                ttemp = tnew - (1.0 - frac)
                ### take integer part of ttemp and sum up
#                aint = sum(counts[before_ind+1:before_ind+1+int(math.floor(ttemp))])
                ### fracind is the index of the last old bin that is split up between the current new bin and the next
#                fracind = int(before_ind + 1 + math.floor(ttemp))
                # print "fracind 2 : " + str(fracind)
                ### redefine frac
#                frac = ttemp - math.floor(ttemp)
                #print "frac: " + str(frac)
#                if frac < 1.0e-10:
#                    frac = 0
                ### if frac is not zero, calculate the part of the old bin that will be in the current new bin
 #               if not frac == 0:
 #                   afrac2 = frac*counts[fracind]
                    #print "afrac2: " + str(afrac2)
#                    cbin.append(afrac1 + aint + afrac2)
#                else:
#                    cbin.append(afrac1+aint)
                #print "cbin: " + str(cbin[-1])
#                before_ind = fracind
                #print "before_ind 2:" + str(before_ind)

#    if method in ['mean', 'avg', 'average', 'arithmetic mean']:
#        print "Method " + method + " chosen."
#        cbin = [c/tnew for c in cbin]
#    elif method not in ['sum']:
#        raise Exception("Method for summing or averaging not recognized. Please enter either 'sum' or 'mean'.")

#    return bintime, cbin, dt

