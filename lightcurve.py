#!/usr/bin/env python
#####################
#
# Class definition for the light curve class. 
# Used to create light curves out of photon counting data
# or to save existing light curves in a class that's easy to use.
#
#

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('PS')
# import matplotlib.pyplot as plt

from generaltools import ask_limits, ask_input, conversion
import sys
import numpy
import math
import numpy as np
import fractions
import scipy.optimize




### some functions to fit lightcurve profiles
def gauss(x, mu, beta, norm):
    return [norm*1.0/(beta * np.sqrt(2.0 * np.pi))* np.exp(-(y - mu)**2.0 / (2.0 * beta**2.0)) for y in x]


def sin(x, a,b,c,d):
    return [b*np.sin(a*y-c)+d for y in x ]



#### BRAND-NEW CLASS IMPLEMENTATION!!!!

class Lightcurve(object):
    def __init__(self, time, counts = None, timestep=1.0, tseg=None, verbose = False, tstart = None):

        if counts == None:
            if verbose == True:
                print "You put in time of arrivals."
                print "Time resolution of light curve: " + str(timestep)
            ### TOA has a list of photon times of arrival
            self.toa = time
            self.ncounts = len(self.toa)
            self.tstart = tstart
            self.makeLightcurve(timestep, tseg = tseg,verbose=verbose)
            
        else:
            self.time = np.array(time)
            self.counts = np.array(counts)
            self.res = time[1] - time[0]
            self.countrate = [t/self.res for t in self.counts]
            self.tseg = self.time[-1] - self.time[0] + self.res

    def makeLightcurve(self, timestep, tseg=None, verbose=False):

        ### if self.counts exists, this is already a light curve, so abort
        try:
            self.counts
            raise Exception("You can't make a light curve out of a light curve! Use rebinLightcurve for rebinning.")
        except AttributeError:

            ## tstart is an optional parameter to set a starting time for the light curve
            ## in case this does not coincide with the first photon
            if self.tstart == None:
                ## if tstart is not set, assume light curve starts with first photon
                tstart = self.toa[0]
            else:
                tstart = self.tstart
            ### number of bins in light curve

            ## compute the number of bins in the light curve
            ## for cases where tseg/timestep are not integer, computer one
            ## last time bin more that we have to subtract in the end
            if tseg:
                timebin = np.ceil(tseg/timestep)
                frac = (tseg/timestep) - int(timebin - 1)
            else:
                timebin = np.ceil((self.toa[-1] - self.toa[0])/timestep)
                frac = (self.toa[-1] - self.toa[0])/timestep - int(timebin - 1)
            print('tstart: ' + str(tstart))

            tend = tstart + timebin*timestep

            ### make histogram
            ## if there are no counts in the light curve, make empty bins
            if self.ncounts == 0:
                print("No counts in light curve!")
                timebins = np.arange(timebin+1)*timestep + tstart
                counts = np.zeros(len(timebins)-1)
                histbins = timebins
                self.res = timebins[1] - timebins[0]
            else:
                timebins = np.arange(timebin+1)*timestep + tstart
                counts, histbins = np.histogram(self.toa, bins=timebin, range = [tstart, tend])
                self.res = histbins[1] - histbins[0]

            #print("len timebins: " + str(len(timebins)))
            if frac > 0.0:
                self.counts = np.array(counts[:-1])
            else:
                self.counts = np.array(counts) 
            ### time resolution of light curve
            if verbose == True:
                print "Please note: "
                print "You specified the time resolution as: " + str(timestep)+ "."
                print "The actual time resolution of the light curve is: " + str(self.res) +"."

            self.countrate = self.counts/self.res
            self.time = np.array([histbins[0] + 0.5*self.res + n*self.res for n in range(int(timebin))])
            if frac > 0.0:
                self.time = np.array(self.time[:-1])
            else:
                self.time = self.time
            self.tseg = self.time[-1] - self.time[0] + self.res

    def saveLightcurve(self, filename):
        """ This method saves a light curve to file. """
        lfile = open(filename, 'w')
        lfile.write("# time \t counts \t countrate \n")
        for t,c,cr in zip(self.time, self.counts, self.countrate):
            lfile.write(str(t) + "\t" + str(c) + "\t" + str(cr) + "\n")
        lfile.close()

    def plot(self, filename, plottype='counts'):
        if plottype in ['counts']:
            plt.plot(self.time, self.counts, lw=3, color='navy', linestyle='steps-mid')
            plt.ylabel('counts', fontsize=18)
        elif plottype in ['countrate']:
            plt.plot(self.time, self.countrate)
            plt.ylabel('countrate', fontsize=18)
        plt.xlabel('time [s]', fontsize=18)
        plt.title('Light curve for observation ' + filename)
        plt.savefig(str(filename) + '.ps')
        plt.close()

    def rebinLightcurve(self, newres, method='sum', verbose = False, implementation="new"):
        ### calculate number of bins in new light curve
        nbins = math.floor(self.tseg/newres)+1
        self.binres = self.tseg/nbins
        print "New time resolution is: " + str(self.binres)

        if implementation in ["o", "old"]:
            self.bintime, self.bincounts, self.binres = self._rebin(self.time, self.counts, nbins, method, verbose=verbose)
        else:
            print("I am here")
            self.bintime, self.bincounts, self.binres = self._rebin_new(self.time, self.counts, newres, method)

    def bkgestimate(self, tseg, loc='both'):
       
        tmin = np.array(self.time).searchsorted(self.time[0]+tseg)
        tmax = np.array(self.time).searchsorted(self.time[-1]-tseg)
        cmin = np.mean(self.counts[:tmin])
        cmax = np.mean(self.counts[tmax:])
        if loc == 'both':
            print("The mean counts/bin before the burst is: " + str(cmin))
            print("The mean counts/bin after the burst is: " + str(cmax))
            print("The combined mean counts/bin is : " + str(np.mean([cmin, cmax])))
            self.meanbkg = np.mean([cmin, cmax])
        elif loc == 'before':
            print("The mean counts/bin before the burst is: " + str(cmin))
            self.meanbkg = cmin
        elif loc == 'after':
            print("The mean counts/bin after the burst is: " + str(cmax))
            self.meanbkg = cmax
        return


    def removebkg(self, tseg, loc='both'):
        self.bkgestimate(tseg, loc=loc)
        counts = self.counts - self.meanbkg
        zeroinds = np.where(counts <= 0.0)[0]
        time = np.array([t for i,t in enumerate(self.time) if not i in zeroinds ])
        counts = np.array([c for i,c in enumerate(counts) if not i in zeroinds ])

        self.ctime = time
        self.ccounts = counts
        return

    ### add Poisson noise to a light curve
    ### this is of some use for artificially created light curves
    def addpoisson(self):
        pcounts = np.array([np.random.poisson for x in self.ctype])
        pcountrate = pcounts/self.res
        self.counts = pcounts
        self.countrate = pcountrate


    ### chop up light curve in pieces and save each piece in
    ### a separate light curve object
    ## len [float]: length of segment (in seconds)
    ## overlap [float, < 1.0]: overlap between segments, in seconds
    def moving_bins(self, timestep=1.0, length=1.0, overlap=0.1):

        #print('self.toa' + str(len(self.toa)))
        ### number of light curves
        nbins = int(math.floor((self.tseg-2.0*overlap)/length))
        print("<<< nbins: " + str(nbins)) 
        try: 
            tstart = self.toa[0]
        except AttributeError:
            raise Exception('No time of arrivals given! Cannot chop up light curve!')

        lcs = []
        tend = 0.0

        while tend <= self.toa[-1] :
            tend = tstart + length 
            stind = self.toa.searchsorted(tstart)
            #print("<<<--- start index : " + str(stind))
            eind = self.toa.searchsorted(tend)
            #print("<<<--- end index: " + str(eind))
            tnew = self.toa[stind:eind]
            #print("<<<--- self.toa: " + str(self.toa[-1]))
            #print("<<<--- tend: " + str(tend))
            if len(tnew) == 0:
                if self.toa[-1] - tend > 0.0:
                    print("tend smaller than end of light curve. Continuing ...")
                    tstart = tend - overlap
                    continue
                else:
                    break
            lcs.append(Lightcurve(tnew, timestep=timestep, tseg=length)) 
            tstart = tend - overlap

        return lcs

    def fitprofile(self, func, p0=None):
        if not p0:
            
            p0 = [10.0, 0.01, 0.01]

        #varobs = np.sum(self.countrate)
        #varmod = np.sum(func(self.time, *p0))
        #renorm = varobs/varmod
        #p0[1] = p0[1] + renorm

        popt, pcov = scipy.optimize.curve_fit(func, self.time, self.counts, p0=p0, maxfev = 50000)
        stderr = np.sqrt(np.diag(pcov))

        print("The best-fit parameters for the FRED model are: \n")
        #print("normalization A = \t" + str(popt[0]) + " \t +/- " + str(stderr[0]))
        #print("burst rise tau1 = \t" + str(popt[1]) + " \t +/- " + str(stderr[1]))
        #print("burst decay tau2 = \t" + str(popt[2]) + " \t +/- " + str(stderr[2]))

        bestfit = func(self.time, *popt)
        newfit = np.where(np.log(bestfit) > 100.0, 1.0e-100, bestfit)
        fitparams = {"popt":popt, "cov":pcov, "err":stderr, "mfit":newfit}

        return fitparams

    def _rebin_new(self, time, counts, dtnew, method='sum'):


        try:
            step_size = float(dtnew)/float(self.res)
        except AttributeError:
            step_size = float(dtnew)/float(self.df)

        output = []
        for i in numpy.arange(0, len(counts), step_size):
            total = 0
            #print "Bin is " + str(i)

            prev_frac = int(i+1) - i
            prev_bin = int(i)
            #print "Fractional part of bin %d is %f"  %(prev_bin, prev_frac)
            total += prev_frac * counts[prev_bin]

            if i + step_size < len(time):
                # Fractional part of next bin:
                next_frac = i+step_size - int(i+step_size)
                next_bin = int(i+step_size)
                #print "Fractional part of bin %d is %f"  %(next_bin, next_frac)
                total += next_frac * counts[next_bin]

            #print "Fully included bins: %d to %d" % (int(i+1), int(i+step_size)-1)
            total += sum(counts[int(i+1):int(i+step_size)])
            output.append(total)

        tnew = np.arange(len(output))*dtnew + time[0]
        if method in ['mean', 'avg', 'average', 'arithmetic mean']:
            cbinnew = output
            cbin = np.array(cbinnew)/float(step_size)
        elif method not in ['sum']:
            raise Exception("Method for summing or averaging not recognized. Please enter either 'sum' or 'mean'.")
        else:
            cbin = output

        return tnew, cbin, dtnew


    ### this method rebins a light curve to a new number of bins 'newbins'
    def _rebin(self, time, counts, newbins, method = 'sum', verbose = False):

        ### nbins is the number of bins in the new light curve
        nbins = int(newbins)
        ### told is the _old_ time resolution
        told = time[1] - time[0]

        ### tseg: length of the entire segment
        tseg = time[-1] - time[0] #+ told
        #print "tseg: " + str(tseg)

        if verbose == True:
            print "nbins: " + str(nbins)
            print "told: " + str(told)
            print "tseg: " + str(tseg)

        ### move times to _beginning_ of each bin
        btime = np.array(time) - told/2.0

        ### dt: new time resolution
        dt = float(tseg)/float(nbins)

        ### check whether old time resolution is larger than new time resolution
        if dt <= told:
            if verbose == True:
                print "Old time resolution bigger than new time resolution."
                print "That's not implemented yet. Returning power spectrum with original resolution."
            return time, counts, told


        ### tnew is the ratio of new to old bins
        tnew = dt/told

        #print "dt: " + str(dt)
        #print "told: " + str(told)
        #print "tnew: " + str(tnew)

        ### new array with bin midtimes
        bintime = [time[0] + 0.5*dt + t*dt for t in range(nbins)]

        ### this fraction is useful, because I can define a number of old bins until the 
        ### boundaries of the old and new bins match up again
        ### this makes it easier in some cases
        tnewfrac = fractions.Fraction(tnew)

        top = tnewfrac.numerator
        bottom = tnewfrac.denominator
        
        #print "top: " + str(top)
        #print "bottom: " + str(bottom)


        ### if the fraction turns out insanely big (happens due to rounding errors), then I do
        ### only one iteration (since tseg/tseg = 1)
        if top > tseg:
            top = tseg
            bottom = nbins
#            print "top: " + str(top)
#            print "bottom: " + str(bottom)
        cbin = []

        ### now iterate over all cycles
#        print "int(tseg/top): " + str(int(nbins/bottom))
#        print("nbins: " + str(nbins)) 

        for i in range(int(nbins/bottom)):

        ### I need this index to remember where I left off during the iteration
            before_ind = 0
#            print "i: " + str(i)
            ### for each cycle iterate through the number of new bins in that cycle
            for j in range(bottom):
                # print "j: " + str(j)
                ### in the first round, start at the lower edge of the bin:
                if before_ind == 0:
                    #print "tnew: " + str(tnew)
                    ## this is the first index to use
                    i0 = int(i*top)
                    #print "i0: " + str(i0)
                    ### first I sum up all complete old bins in that new bin
                    aint = sum(counts[i0:int(i0+math.floor(tnew))])
                    #print "lower index: " + str(i0)
                    #print "upper index: " + str(int(i0+math.floor(tnew)))
                    #print "values to sum: " + str(counts[i0:int(i0+math.floor(tnew))])

                    ### then I set the index of the old bin that is _not_ completely in the new bin
                    fracind = int(i0 + math.floor(tnew) )
                    #print "fracind 1 : " + str(fracind)


                    ### frac is the fraction of the old bin that's in the new bin
                    frac = tnew - math.floor(tnew)
                    #print "tnew fractional part: "  + str(tnew- math.floor(tnew))

                    ### if frac is not zero, then compute fraction of counts that goes into my new bin
                    if frac < 1.0e-10:
                        frac =0
                    if not frac == 0:
                        afrac = frac*counts[fracind]
                        #print "afrac: " + str(afrac)
                        cbin.append(aint + afrac) ### append to array with new counts
                    else:
                        cbin.append(aint)
                        #print "cbin: " + str(cbin[-1])

                    ### reset before_ind for next iteration in j
                    before_ind = fracind
                    #print "before_ind 1 : " + str(before_ind)
                else:

                    ### This new bin doesn't start at the same position as the old bin, hence I start with the fraction
                    ### afrac1 is the rest of the preceding old bin that was split up
                    afrac1 = (1.0 - frac)*counts[before_ind]
                    # print "afrac1: " + str(afrac1)
                    ### 1.0-frac of the bin already done, so define new length for the rest: ttemp 
                    ttemp = tnew - (1.0 - frac)
                    ### take integer part of ttemp and sum up
                    aint = sum(counts[before_ind+1:before_ind+1+int(math.floor(ttemp))])
                    ### fracind is the index of the last old bin that is split up between the current new bin and the next
                    fracind = np.int(before_ind + 1 + math.floor(ttemp))
                    #print "fracind 2 : " + str(fracind)
                    ### redefine frac
                    frac = ttemp - math.floor(ttemp)
                    #print "frac: " + str(frac)
                    if frac < 1.0e-10:
                        frac = 0
                    ### if frac is not zero, calculate the part of the old bin that will be in the current new bin
                    if not frac == 0:
                        #print("fracind2: " + str(fracind))
                        afrac2 = frac*counts[int(fracind)]
                        #print "afrac2: " + str(afrac2)
                        cbin.append(afrac1 + aint + afrac2)
                    else:
                        cbin.append(afrac1+aint)
                    #print "cbin: " + str(cbin[-1])
                    before_ind = fracind

        if method in ['mean', 'avg', 'average', 'arithmetic mean']:
            cbinnew = cbin
            cbin = [c/tnew for c in cbinnew]
        elif method not in ['sum']:
            raise Exception("Method for summing or averaging not recognized. Please enter either 'sum' or 'mean'.")
        return bintime, cbin, dt

###############################################################

###############################################################
##### FUNCTIONS ###############################################
###############################################################


#### ADD NOISE TO A LIGHT CURVE ##############################
#
#
#
#
#
#
def addnoise(lc):
   
    #time = lc.time
    #counts = lc.counts

    for i,t in enumerate(lc.time):
        pmean = lc.counts[i] 
        lc.counts[i] = np.random.poisson(pmean)

    return




# MAKE LIGHT CURVES FROM TTE DATA IN DIFFERENT ENERGY RANGES ######################
#
# !!! 21/03/11: made changes from taking a number of bins to specifying only emin and emax bounds
# NEED TO IMPLEMENT THAT CHANGE PROPERLY!!!
#
#
#
#
# This script takes energy ranges specified by the user, separates the tte data into different
# lists and makes a light curve from each.
#
# REQUIRES: numpy, math
#
# INPUT: 1) burst number (bn)
#	 2) detector number (detec)
#	 3) number of energy bins
#	 [4) time resolution of the light curve (in seconds); currently queried for in main()]
#
# NOTE: this script uses the output of channeltoenergy.py --> won't work on channel data!
#
# OUTPUT: 
#
#
#

#def energybins(bins, tnew, evt_en, emin, emax):
#    ebin, en, eb=[], [], []
#    print "This is the minimum energy chosen: " + str(bins[0])
#    print "And this is the minimum energy emin[4]: " + str(emin[4])
#    if float(bins[0]) < float(emin[4]):
#        print "The minimum energy specified is smaller than the low-energy cut-off for reliable data. Setting minimum energy to " + str(emin[4])
#        bins[0] = emin[4]
#    if float(bins[1]) > float(emax[-2]):
#        print "The maximum energy specified is larger than the high-energy cut-off for reliable data. Setting maximum energy to " + str(emax[-2])
#        bins[1] = emax[-2]
#    ttetemp, evttemp = [], []
#    for j,temp in enumerate(tnew):
#        if bins[0] <= evt_en[j] < bins[1]:
#            ttetemp.append(temp)
#            evttemp.append(evt_en[j])
#        else: continue
#    print "This is a test 2"
#    tte = {'tebin': ttetemp, 'evtebin':evttemp}
#    print "tte keys: "  + str(tte.keys)
#    return tte
 


# THIS CODE MAKES A NUMBER OF ENERGY BINS, FOR NOW I ONLY WANT ONE WHERE I CAN SPECIFY THE MIN AND MAX ENERGY
#    te=(emax[-2] - emin[4])/(int(bins)+1)  # leave first four and last channel free
#    for i in range(int(bins)+1):
#        if i==0: 
#            enow=emin[4]
#        else:
#            enow=enow+te
#        ebin.append(enow)
#    print "This is emin: " + str(emin[4]) + " and this emax: " + str(emax[-2])
#    print ebin
#    return ebin

#def lc(tnew, timestep):
#    timebin=math.floor((tnew[-1]-tnew[0])/timestep)+1
#    counts,histbins=numpy.histogram(tnew, bins=timebin)
#    timestepnew = histbins[1] - histbins[0]
#    lctimes = [histbins[0]+((0.5+n)*timestep) for n in range(int(timebin))]
#    lcdict = {'timestep': timestepnew, 'timebin':timebin, 'histbins':histbins, 'lctimes':lctimes, 'counts': counts}
#    return lcdict

#def main():
#    filename=sys.argv[1]
#    bn=sys.argv[1]
#    print "You chose burst number " + str(bn)
#    detec=sys.argv[2]
#    print "You chose detector number: " + str(detec)
#    filename='bn'+ str(bn) + '_n' + str(detec) + '_tte_energy.dat'
#    print "The input file is: " + filename
#    if sys.argv[3] == 0:
#        bins=1
#    else: bins=[sys.argv[3], sys.argv[4]]
#    print "The number of energy bins is: " + str(bins)
#    ttelist= conversion(filename)
#    tnew=ttelist[0]
#    events=ttelist[1]
#    channelconv = conversion('channelconv.dat')
#    emin=channelconv[1]
#    emax=channelconv[2]
#    emid=channelconv[3]
#    ebin = energybins(bins, emin, emax)
#    tte,lcarray = lc(tnew, events, ebin, bins)
#    for i, my_lc in enumerate(lcarray):
#        ttefile=open(str(bn) + '_n' + str(detec) + '_tte_' + str(int(ebin[i])) + 'keVto' + str(int(ebin[i+1])) + 'keV.dat', 'w')
#        ttefile.write('#[time] \t [event]')
#        ttefile.write('\n')
#        for tbin, evtbin in zip(tte[i][0], tte[i][1]):
#            ttefile.write(str(tbin) + "\t" + str(evtbin))
#            ttefile.write('\n')
#        myfile=open(str(bn) + '_n' + str(detec) + '_lc_' + str(int(ebin[i])) + 'keVto' + str(int(ebin[i+1])) + 'keV.dat', 'w')
#        myfile.write('#[time bin] \t [events in energy]')
#        myfile.write('\n')
#        for histbin, n in zip(my_lc[0], my_lc[1]):
#            myfile.write(str(histbin) + "\t" + str(n))
#            myfile.write('\n')
#        myfile.close()
#        ttefile.close()
#    return

#if __name__ == "__main__":
#    main()

