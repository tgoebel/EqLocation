#!/usr/bin/python3
"""
    - run 2_cd11_2_mseedDy.py first
    - this script loads resulting output and creates basci plots
"""
#=================================================================================
#                       modules
#=================================================================================
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
from obspy.core import read
from obspy.core import UTCDateTime

#=================================================================================
#                       files and parameters
#=================================================================================
mseed_dir= '/media/tgoebel/Data/seis/BlueMountain/mseed'
l_Sta    = ['BM01', 'BM02', 'BM03', 'BM04', 'BM06', 'BM07', 'BM09', 'BM10']
l_Sta    = ['BM03']# 'BM02',
# bandpass filter
f_min, f_max = .1, 50
test_plot = True
#tmin     = UTCDateTime('2019-05-02')
tmin     = UTCDateTime('2019-06-25')
tmax     = UTCDateTime('2019-06-26')

for sta in l_Sta:
    print(  '--------------------- station', sta, '----------------------------')
    curr_t    = tmin
    n_noReads = 0
    #===============================2=================================================
    #                      load mseed files
    #=================================================================================
    while (tmax-curr_t) > 0:
        t_str = '%04d%02d%02d'%(curr_t.year,curr_t.month,curr_t.day)
        #write to mseed
        mseed_file = '%s_%s.mseed'%( sta, t_str)
        SEIS = read( '%s/%s/%s'%( mseed_dir, sta, mseed_file))
        ###B### detrend and filtering, bandpass
        SEIS.detrend( type = 'linear')# linear or simple or demean
        #SEIS.filter( 'bandpass', freqmin=f_min, freqmax=f_max, corners = 2)
        print( SEIS[0].stats)
        #SEIS[0].plot()
        if test_plot == True:
            a_tr1 = SEIS[0]
            dt = a_tr1.stats.delta
            a_t_hr = np.arange( 0, dt*a_tr1.stats.npts, dt)#/3600.
            n_axis = int( a_t_hr[-1])+1
            # plot in hourly axis increments
            #plt.figure(1, figsize=(6,14))
            for i_ax in range( n_axis):
                #ax = plt.subplot( n_axis, 1, i_ax+1)#, sharex = ax)
                plt.figure()
                ax = plt.subplot( 111)
                i1,i2 = int(i_ax*(3600./dt)), int((i_ax+1)*3600./dt)
                ax.plot( a_t_hr[i1:i2], a_tr1.data[i1:i2])
                ax.set_xlabel( tmin)
                plt.show()
        curr_t += 3600*24

