#!/usr/bin/env python

# INTERPOLATE JPL EPHEMERIS POSITIONS
#
# Fairly straightforward translation from IDL into python
# of Craig Markwardt's IDL routine jplephinterp.PRO
#
#
#
import numpy as np
from jplephread import jplephread
#from readposhist import poshist
import pyfits
import sys
from math import floor
import time as tsystem

def jplinterp(obj, info, time, tbase, coeffs):
    tstart= tsystem.clock()
    if obj == 3:
        print "Translating from earth-moon barycenter to earth"
        # EARTH-MOON EPHEMERIDES
        xem, yem, zem, vxem, vyem, vzem = jplephcalc(time, info, coeffs, tbase, obj)
        xmo, ymo, zmo, vxmo, vymo, vzmo = jplephcalc(time, info, coeffs, tbase, 10)
        emrat = info['emrat']
        x, y, z, vx, vy, vz = [], [], [], [], [], []
        for i, (xe, xm, ye, ym, ze, zm, vxe, vxm, vye, vym, vze, vzm) in enumerate(zip(xem, xmo, yem, ymo, zem, zmo, vxem, vxmo, vyem, vymo, vzem, vzmo)):
            x.append(xe - emrat*xm)
            y.append(ye - emrat*ym)
            z.append(ze - emrat*zm)     
            vx.append(vxe - emrat*vxm)
            vy.append(vye - emrat*vym)
            vz.append(vze - emrat*vzm)
        tend = tsystem.clock()
        print "JPLINTERP Execution time, EARTH: " + str(tend - tstart)
        return x,y,z,vx,vy,vz
    elif obj == 11:
        print "Calculating solar ephemeris"
        x, y, z, vx, vy, vz = jplephcalc(time, info, coeffs, tbase, 11)
        tend = tsystem.clock()
        print "JPLINTERP Execution time, SUN: " + str(tend - tstart)
        return x,y,z,vx,vy,vz
    else:
        print "Your object has not been implemented yet." 

def jplephcalc(time, info, coeffsnew, tbase, objold):
# time is time of interest, relative to tbase (which units?), actual interpulation time is (t+tbase)
# objold:  sun = 11, earth = 3, solarbary, ssb = 12, earthbary = 13
# I think I only need earth and sun for my analysis
# these numbers are to be consistent with idl and fortran code
# the actual index is objold-1
#    tbase=mjdrefi+ 2400000.5
#    time = []
#    for i in phlist:
#        time.append((i + sec_offset)/8.64e4)
    obj=objold-1 # number of object; change for earth
    nc = info['ncoeff'][obj]
    ns = info['nsub'][obj]
    dt = info['timedel']
    nr = info['jdrows']
    jd0 = info['jdlimits'][0] - tbase
    jd1 = info['jdlimits'][1] - tbase

    if obj==12:
        print "You chose SUNBARY. That's not incorporated yet. Please come back later."
        return
    else:
        ii1 = info['ptr'][obj]-1
        ii2 = ii1 + nc*ns*3.0 - 1
        coeffssmall=[]
        for i in coeffsnew:
            coeffssmall.append(i[ii1:ii2+1])
        coeffs=np.reshape(coeffssmall, (nr, ns, 3, nc))

    #decide which interval and subinterval we're in (WHAT DOES THAT MEAN?)
    # Vectorized version: using numpy array
    tint = (np.array(time)-jd0)/dt
    ieph = np.floor(tint)
    tint = (tint - ieph)*ns
    nseg = np.floor(tint)
    tseg = 2*(tint - nseg) - 1

#    for i,j in enumerate(time):
#        tintnow=(j-jd0)/dt
#        ieph.append(floor(tintnow))
#        tint.append((tintnow-ieph[i])*ns)
#        nseg.append(np.floor(tint[i]))
#        tseg.append(2.0*(tint[i] - nseg[i]) - 1)

    # optimization
    #try:
    #    len(ieph)
    mini=[min(ieph), max(ieph)]
#    print "mini: " + str(mini)
    minn=[min(nseg), max(nseg)]
#    print "minn: " + str(minn)
    if mini[0] == mini[1] and minn[0] == minn[1]:
        ieph=int(ieph[0])
        nseg=int(nseg[0])
        i0 = int(0)
    else:
        ieph = np.array([int(i) for i in ieph])
        nseg = np.array([int(i) for i in nseg])
        i0 = np.array([int(i) for i in np.zeros(len(ieph))]) # for now I assume that ieph is a scalar, but I might have to assume it's not in the future!

    p0 = np.ones(len(tseg))
    p1 = tseg
    v0 = np.zeros(len(p1))
    v1 = np.ones(len(p1))
 
    xtemp, ytemp, ztemp, vxtemp, vytemp, vztemp, p2, p3, v3 = np.zeros(len(p1)), np.zeros(len(p1)), np.zeros(len(p1)), np.zeros(len(p1)), np.zeros(len(p1)), np.zeros(len(p1)), np.zeros(len(p1)), np.zeros(len(p1)), np.zeros(len(p1))

    #except TypeError:
    #    p0 = 1.0
    #    p1 = tseg
    #    v0 = 0.0
    #    v1 = 1.0 
    #    i0 = 0.0
    #    xtemp, ytemp, ztemp, vxtemp, vytemp, vztemp, p2, p3, v3 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    
    v2 = 4.0*tseg
    tt = 2.0*tseg

    #v2, tt = [], []
    #for i in tseg:
    #    v2.append( 4.0*i)
    #    tt.append( 2.0 * i)
       

    i1 = i0 + 1
    i2 = i1 + 1

    vfac = 2.0*ns/dt
    #lenp = np.zeros(len(p1))
    #xtemp, ytemp, ztemp, vxtemp, vytemp, vztemp, p2, p3, v3 = np.zeros(len(p1)), np.zeros(len(p1)), np.zeros(len(p1)), np.zeros(len(p1)), np.zeros(len(p1)), np.zeros(len(p1)), np.zeros(len(p1)), np.zeros(len(p1)), np.zeros(len(p1))
    i=0
    if nc%2 ==0:
        ncnew=nc-1
    else:
        ncnew=nc+1
    for i in range(0,ncnew,2):
        if i == nc-1:
            p1 = np.zeros(len(p1))
            v1 = np.zeros(len(p1))
#        for j in i0:  ### Only need this if size(ieph) > 1
        ii = i0 + i 
        jj = i0 + min(i+1, nc-1)
# Note: xtemp, vxtemp, p2 etc are numpy arrays --> these are vector operations!

#        print "ieph: " + str(ieph)
#        print "nseg: " + str(nseg)
#        print "i0: " + str(i0)
#        print "ii: " + str(ii)
#        print "jj: " + str(jj)
        xtemp=xtemp + coeffs[ieph, nseg, i0, ii]*p0 + coeffs[ieph, nseg, i0, jj]*p1
        ytemp=ytemp + coeffs[ieph, nseg, i1, ii]*p0 + coeffs[ieph, nseg, i1, jj]*p1
        ztemp=ztemp + coeffs[ieph, nseg, i2, ii]*p0 + coeffs[ieph, nseg, i2, jj]*p1
        vxtemp=vxtemp + coeffs[ieph, nseg, i0, ii]*v0 + coeffs[ieph, nseg, i0, jj]*v1
        vytemp=vytemp + coeffs[ieph, nseg, i1, ii]*v0 + coeffs[ieph, nseg, i1, jj]*v1
        vztemp=vztemp + coeffs[ieph, nseg, i2, ii]*v0 + coeffs[ieph, nseg, i2, jj]*v1
        p2 = tt*p1 - p0
        p3= tt*p2 - p1
        v2= tt*v1 - v0 + 2.0*p1
        v3 = tt*v2 - v1 + 2.0*p2
#        for k, pol in enumerate(p1):
#            xnow = coeffs[ieph, nseg, i0, ii]*p0[k] + coeffs[ieph, nseg, i0, jj]*pol #calculate x for each element in p1, v1
#            xtemp[k] = xtemp[k] + xnow  # append to temporary x-array
#            ynow = ytemp[k] + coeffs[ieph, nseg, i1, ii]*p0[k] + coeffs[ieph, nseg, i1, jj]*pol #calculate y for each element in p1, v1
#            ytemp[k] = ynow # append to temporary y-array
#            znow = ztemp[k] + coeffs[ieph, nseg, i2, ii]*p0[k] + coeffs[ieph, nseg, i2, jj]*pol #calculate z for each element in p1, v1
#            ztemp[k] = znow # append to temporary z-array
#            vxnow = vxtemp[k] + coeffs[ieph, nseg, i0, ii]*v0[k] + coeffs[ieph, nseg, i0, jj]*v1[k]
#            vxtemp[k] = vxnow  # velocities in km/s
#            vynow = vytemp[k] + coeffs[ieph, nseg, i1, ii]*v0[k] + coeffs[ieph, nseg, i1, jj]*v1[k]
#            vytemp[k] =vynow
#            vznow = vztemp[k] + coeffs[ieph, nseg, i2, ii]*v0[k] +  coeffs[ieph, nseg, i2, jj]*v1[k]
#            vztemp[k] = vznow
#            p2[k] = tt[k]*pol - p0[k]
#            p3[k]= tt[k]*p2[k] - pol
#            v2[k]= tt[k]*v1[k] - v0[k] + 2.0*pol
#            v3[k] = tt[k]*v2[k] - v1[k] + 2.0*p2[k]      
        p0 = p2.copy()
        p1 = p3.copy()
        v0 = v2.copy()
        v1 = v3.copy()
    vx=vxtemp*vfac/8.64e4
    vy=vytemp*vfac/8.64e4
    vz=vztemp*vfac/8.64e4
    return xtemp, ytemp, ztemp, vx, vy, vz

#def main():
#    date=str(sys.argv[1])
#    dateno=date[8:]+date[3:5]+date[:2]
#    print "You chose position information for date: " + date
#    filename='glg_poshist_all_' + dateno + '_v00.fit'
#    print filename
#    phlist = poshist(filename)
#    info, coeffs = ephread(phlist['jdlimits'])
#    obj = int(sys.argv[2])
#    print "You have chosen object: " + str(obj)
#    time=np.array([])
#    for i in phlist['phtime']:
#    time=(phlist['phtime'] + phlist['sec_offset'])/8.64e4
#        time.append((i+phlist['sec_offset'])/8.64e4)	#time in fraction JD from base_jd
#    print time
#    tbase = phlist['mjdrefi']+2400000.5
#    x, y, z, vx, vy, vz = jplinterp(obj, info, time, tbase, coeffs)
    # For some reason, jplinterp gives out a *list* of a numpy array. Why? Dunno, need to figure out!
#    x=x[0]
#    y=y[0]
#    z=z[0]
#    vx=vx[0]
#    vy=vy[0]
#    vz=vz[0]
#    outfile = sys.argv[3]
#    print "The outfile will be named: " + outfile
#    efile=open(outfile, 'w')
#    efile.write('#x \t y \t z \t vx \t vy \t vz ')
#    efile.write('\n')
#    for (a,b,c,d,e,f) in zip(x, y, z, vx, vy, vz):
#        efile.write(str(a) + '\t' + str(b) + '\t' + str(c) + '\t' + str(d) + '\t' + str(e) + '\t' + str(f))
#        efile.write('\n')
#    efile.close()

#if __name__ == "__main__":
#    main()
