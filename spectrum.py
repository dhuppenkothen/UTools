
import numpy
import scipy.optimize

import generaltools as gt


class Spectrum(object):

    def __init__(self, data = None, freq = None, flux = None, unit='erg'):

        if data == None and not freq == None and not flux == None:
            self.freq = np.array(frequency)
            self.flux = np.array(flux)

        elif not data == None and freq == None and flux == None:
            self.read_spectrum(data)

        else:
            self.freq = None
            self.flux = None

    def read_spectrum(self, filename):

        raw_data = gt.conversion(data)
        self.freq = np.array([float(x) for x in raw_data[0]])
        self.flux = np.array([float(x) for x in raw_data[0]])


    def convert_flux(self, new_unit='jansky'):

        if unit.lower() in ['jansky', 'jy']:
            fnew = self.flux*1.0e23

        return fnew


    ### FITTING FUNCTION func TO DATA ####
    def fit_spectrum(func, init_guess=None, sigma = None, method='chisquare'):
 
        if method.lower() in ['chisquare', 'chi2', 'chi']:
            fpar, fcov = scipy.optimize.curve_fit(func, self.freq, self.flux, p0=init_guess, sigma = sigma)
            print("The best-fit parameters with errors are: \n")
  

