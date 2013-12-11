#!/usr/bin/env python

#### READ IN JPL EPHEMERIS FILE
#
# Straightforward translation of JPLEPHREAD.PRO written in IDL
# by Craig Markwardt
#
#
# 
# TO DO: - do all the different error checks like in IDL program
#
#

import pyfits
from math import floor
import time as tsys
import os




def ephread(jdlimits, ephname = 'JPLEPH.200', ephdir='current'):
    tsysstart = tsys.clock()
    # READ IN EPHEMERIS FILE

    ### specify directory of ephemeris file with ephdir keyword
    ### if none specified, assume ephemeris file is in current working directory
    if ephdir == 'current':
        ephdir = os.path.dirname(__file__)
    ### standard ephemeris file used is JPLEPH.200, but other specifications possbile
    ### (not tested, so not sure other files work!)
    ephfile = os.path.join(ephdir, ephname)

    hdulist=pyfits.open(ephfile)
    # READ OUT FIRST EXTENSION AFTER PRIMARY HEADER
    ext1=hdulist[1].data
    cname=ext1.field(0)
    cval=ext1.field(1)
    # DEFINE SOME CONSTANTS
    denum = round(cval[cname == 'DENUM'])
    clightold = float(cval[cname=='CLIGHT'])
    emrat = 1.0/(1.0 + float(cval[cname == 'EMRAT']))
    au = float(cval[cname == 'AU'])
    msolold = float(cval[cname == 'GMS'])
    radsolold= float(cval[cname == 'RADS'])

    # CALCULATE constants in funky units
    x = au/clightold 
    msol = msolold*(x**3.0)/(86400.0**2.0) #G*M_sun in light seconds
    radsol = radsolold/clightold
    clight = clightold * 1000.0  # speed of light in m/s
    # FIRST EXTENSION DONE, NOW SECOND EXTENSION
    ext2=hdulist[2].data
    ephobj = ext2.field(0)
    ephptr = ext2.field(1)
    ephncoeff = ext2.field(2)
    ephnsub = ext2.field(3)
    #THIRD EXTENSION
    ext3=hdulist[3].data
    ext3hdr=hdulist[3].header
    nrows=ext3hdr['NAXIS2']
    tstart=ext3hdr['TSTART']
    tstop=ext3hdr['TSTOP']
    timedel=ext3hdr['TIMEDEL']
    rowlimits=[]
    expand=[-2,2]
    for i in range(2):
         rowlimits.append(floor((jdlimits[i] - tstart)/timedel) + expand[i])


    # READ IN RAW DATA
    coeffs=ext3['ChebCoeffs']
    cols=hdulist[3].columns
    dims=len(coeffs[0])
    coeffsnew=coeffs[rowlimits[0]-1 : rowlimits[1]]
    nr = rowlimits[1] - rowlimits[0]+1
    jdlimits1 = [0,0]
    jdlimits1[0] = (rowlimits[0]-1)*timedel + tstart
    jdlimits1[1] = rowlimits[1]*timedel + tstart
    info = {'nrows': nrows, 'tstart' : tstart, 'tstop': tstop, 'timedel': timedel, 'format': 'BINEPH2FITS', 'denum': denum, 'c': clight, 'emrat': emrat, 'au':au*1000/clight, 'msol': msol, 'sunrad': radsol, 'jdlimits': jdlimits1, 'jdrows': nr, 'objname': ephobj, 'ptr': ephptr, 'ncoeff':ephncoeff, 'nsub':ephnsub}
    tend = tsys.clock()
    print "EPHREAD Execution Time: " + str(tend - tsysstart)
    return info, coeffsnew

# Need to finish this script, return some more values --> read about classes!
 

