UTools: Bits and pieces of code I use a lot

UTools contains a selection of tools for various purposes.# This is everything I need across projects or did not fit elsewhere.

See descriptions below.

Contents:

generaltools: a collection of functions of general usefulness, for example
	- conversion: read in ascii data in columns into a python dictionary
	- getpickle: read out pickled files


lightcurve.py: create light curves from photon arrival times, rebin, plot and manipulate
 
mle.py: Maximum likelihood fitting for periodograms, light curves and, more recently, Gaussian Processes

posterior.py: defines a Posterior object, with subclasses for different types (periodograms, light curves, Gaussian Processes)

priors.py: defines functional forms to be used for priors

simpower.py: make red noise light curves from power spectra using Timmer+Koenig 1995

spectrum.py [UNTESTED]: make energy spectra from time tagged event data with attached energies

ttrig_to_utc.py: convert time stamps in time since trigger (Fermi/GBM data) to UTC

xcor.py [MOSTLY UNTESTED]: make convolutions, cross correlations and cross spectra
 

