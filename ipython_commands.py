import matplotlib.pyplot as plt
from pylab import *

import numpy as np

import generaltools as gt
import lightcurve
import powerspectrum
import mle
import bayes

def read_data(filename, bst, fnyquist = 4096.0, norm='variance'):
    dt = 0.5/fnyquist
    procdata = gt.getpickle(filename)
    evt = procdata['combined']
    evt.filterburst(bst)
    time = np.array([x.time for x in evt.burstphot])
    lc = lightcurve.Lightcurve(time, timestep=dt)
    ps = powerspectrum.PowerSpectrum(lc, norm=norm)

    return time, lc, ps

def fit_data(ps, func=mle.pl, ain=[2.0, 4.0, -5.0], fitmethod='powell', bounds=None, qpo = False, plotname='lrttest'):
    fitspec = mle.MaxLikelihood(ps.freq[1:], ps.ps[1:], fitmethod=fitmethod, obs=True)
    fitparams = fitspec.mlest(func, ain, fitmethod=fitmethod, bounds=bounds)
    if qpo:
        lrt, optpars, qpopars = fitspec.find_qpo(func, fitparams['popt'], plot=True, plotname=namestr)
        return lrt, optpars, qpopars
    else:
        return fitparams


def model_selection(ps, namestr='test', nsim=1000, choose_mod = True, find_qpo = True, bestmod=mle.pl, ain=[2.0, 5.0, -5.0]):
    btest = bayes.Bayes(ps, namestr=namestr)
    if choose_mod:
        ain_bpl = [1.0, ain[1], ain[0], 1.0, ain[2]]
        psfit, fakeper, p_lrt = btest.choose_noise_model(mle.pl, ain, mle.bpl, ain_bpl, nsim=nsim)
        if p_lrt < 0.05:
            bestmod = mle.bpl
            ain = psfit.bplfit['popt']
        else:
            bestmod = mle.pl
            ain = psfit.plfit['popt']
    if find_qpo:
        result = btest.find_periodicity(bestmod, ain, nsim=nsim)
        btest.find_qpo(bestmod, ain, nsim=nsim, plot=True, plotstr=namestr)
    return 

def make_panel_plot(bids, bsts, plotname):


    figprops = dict(figsize=(18.,10.), dpi=128)
    adjustprops = dict(left=0.1, bottom=0.1, right=0.98, top=0.93, wspace=0.2, hspace=0.2)


    f = plt.figure(**figprops)
    ax = f.add_subplot(111)
    ax.spines['top'].set_color("none")
    ax.spines['bottom'].set_color("none")
    ax.spines['left'].set_color("none")
    ax.spines['right'].set_color("none")
    ax.tick_params(labelcolor="w", top="off", bottom="off", left="off", right="off")


    f.subplots_adjust(**adjustprops)

    colnum = 3
    rownum = int(len(bids)/3)
    print("rownum: " + str(rownum))
    if not colnum*rownum == len(bids):
       rownum = rownum + 1

    plotcount = 0


    for x,b in zip(bids, bsts):
        plotcount = plotcount+1

        filename = "tte_bn" + str(x) + "_procdata.dat"
        time, lc, ps = read_data(filename, b)
        lc = lightcurve.Lightcurve(time, timestep=0.005)

        #if plotcount == 1:
        #    ax = subplot(rownum, colnum, plotcount)
        #else:   
        #    subplot(rownum, colnum, plotcount, sharex=ax, sharey=ax)
        ax1 = f.add_subplot(rownum, colnum, plotcount)

        plot(lc.time, lc.countrate, lw=3, c='navy', linestyle="steps-mid")
        plt.axis([min(lc.time), max(lc.time), min(lc.countrate), max(lc.countrate)+np.mean(lc.countrate)])

    ax.set_ylabel("Countrate [counts/s]")
    ax.set_xlabel("time since trigger [s]")
        

    return 


