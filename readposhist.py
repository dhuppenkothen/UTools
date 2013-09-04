#!/usr/bin/env python

## READ POSITION HISTORY FILE FOR GBM ###################################
#
#
# REQUIRES: pyfits, numpy
#
# INPUT: date in formate dd/mm/yy
#
#
#
#

import pyfits
import numpy
import sys
import time as tsys

def poshist(filename):
    tstart = tsys.clock()
    hdulist=pyfits.open(filename)
    mjdrefi=hdulist[0].header['MJDREFI']
    mjdreff=hdulist[0].header['MJDREFF']
    sec_offset = mjdreff*8.64e4
    data=hdulist[1].data
    phtime=data.field(0)	#contains sclk_utc = mission elapsed time in seconds
    position=[data.field('pos_x'), data.field('pos_y'), data.field('pos_z')]
    jdlimits=[]
    for i in [-1, 2]:
        jdlimits.append(phtime[0]/8.64e4 + i + mjdreff + mjdrefi + 2400000.5)
    tbase = mjdrefi+2400000.5
    time=(phtime + sec_offset)/8.64e4
    phlist = {'phtime':phtime, 'time':time, 'basejd':tbase, 'mjdrefi':mjdrefi, 'mjdreff':mjdreff, 'sec_offset':sec_offset, 'position':position, 'jdlimits':jdlimits}
    tend = tsys.clock()
    print "POSHIST Execution Time: " + str(tend-tstart)
    return phlist
 

def main():
    date=str(sys.argv[1])
    dateno=date[8:]+date[3:5]+date[:2]
    print "You chose position information for date: " + date
    filename='glg_poshist_all_' + dateno + '_v00.fit'
    print filename
    mjdrefi, mjdreff, sec_offset, phlist, nph, jdlimits = poshist(filename)
    print str(mjdrefi) + " " + str(mjdreff) 
    print "jdlimits: " + str(jdlimits)
    return


if __name__ == "__main__":
    main()



