import numpy as np
import argparse
 
## convert trigtime to UTC
def ttrig_to_utc(t, trigtime):

    ###trigger time in seconds since midnight that day
    ## step 1: convert MET into days
    ttrig_day = trigtime/(86400.0)

    ## step 2: take fraction of day
    ttrig_frac = ttrig_day - np.floor(ttrig_day)

    ## step 3: convert back into seconds:
    ttrig_sec = ttrig_frac*86400.0

    tsec = t + ttrig_sec - 2.0
    htime = tsec/3600.0
    hours = int(np.floor(htime))
    hfrac  = htime - hours

    mtime = hfrac*60.0 ## = hfrac*3600.0/60.0
    minutes = int(np.floor(mtime))
    mfrac = mtime - minutes

    seconds = mfrac*60.0

    if hours < 10:
       hours = '0' + str(hours)
    if minutes < 10:
       minutes = '0' + str(minutes)

    if np.floor(seconds) < 10:
       second s= '0' + str(seconds)

    utc = str(hours) + ':' + str(minutes) + ':' + str(seconds)

    print('The time in UTC is ' + utc)

    return utc


def main():

    utc = ttrig_to_utc(t, trigtime)
    return


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Converting from time since trigger to UTC.")

    parser.add_argument('-s', '--starttime', dest='t', help='Start time in seconds since trigger.')
    parser.add_argument('-t', '--trigtime', dest='trigtime', help='Time of GBM Observation Trigger.')

    clargs = parser.parse_args()
    t = float(clargs.t)
    trigtime = float(clargs.trigtime)
    main()
    
    


