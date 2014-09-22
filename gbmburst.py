import numpy as np
import cPickle as pickle

import generaltools as gt
import lightcurve
import powerspectrum
import mle
import bayes
from burst import Burst


#### BURST SUBCLASS FOR GBM DATA ###################
#
# Subclass of class Burst to deal with GBM data
# and GBM specific issues
#
#
class GBMBurst(Burst, object):


    def __init__(self, bid, bstart, blength,
                 energies=None,
                 photons=None,
                 events = None,
                 filename=None,
                 instrument="gbm",
                 fnyquist = 4096.0,
                 norm='leahy',
                 fluence = None,
                 epeak = None,
                 ttrig = None):

        ### set burst ID
        self.bid = bid

        ### if photons and filename aren't given, then data comes from procdata file
        if photons ==None and not filename:
           filename = "tte_bn" + str(bid) + "_procdata.dat"
       
        Burst.__init__(self, bstart, blength,
                 energies,
                 photons,
                 events,
                 filename,
                 instrument,
                 fnyquist,
                 norm,
                 fluence = fluence,
                 epeak = epeak,
                 ttrig = ttrig)
        
        return


        def read_data(self, filename, filetype='ascii', det="combined"):

            Burst.read_data(self, filename, type=filetype)

            evt = self.photons[det]
            self.photons = np.array([x.time for x in evt.photons])
            self.energies = np.array([x.energy for x in evt.photons])
            return            
 

 

