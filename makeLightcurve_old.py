


    def makeLightcurve(self, timestep, tseg=None, verbose=False):
        try:
            self.counts
            raise Exception("You can't make a light curve out of a light curve! Use rebinLightcurve for rebinning.")
        except AttributeError:
            ### number of bins in light curve
            if tseg:
                timebin = np.round(tseg/timestep)
            else:
                timebin = np.floor((self.toa[-1] - self.toa[0])/timestep) + 1
            print('timebin: ' + str(timebin))
            ### make histogram
            timebins = np.arange(timebin+1)*timestep + self.toa[0]
            counts, histbins = np.histogram(self.toa, bins=timebins)
            self.counts = np.array(counts)
            ### time resolution of light curve
            self.res = histbins[1] - histbins[0]
            if verbose == True:
                print "Please note: "
                print "You specified the time resolution as: " + str(timestep)+ "."
                print "The actual time resolution of the light curve is: " + str(self.res) +"."

            self.countrate = self.counts/self.res
            self.time = np.array([histbins[0] + 0.5*self.res + n*self.res for n in range(int(timebin))])
            self.tseg = self.time[-1] - self.time[0] + self.res

