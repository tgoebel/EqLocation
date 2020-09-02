#!/usr/bin/python3
"""
    - plot waveforms across all stations, either only
    one channel or all channels
      (use: s_com = 'all, or s_com = 'Z'

    - run 2_cd11_2_mseedDy.py first
    - this script loads resulting output and creates basic plots
"""
#=================================================================================
#                       modules
#=================================================================================
import matplotlib as mpl
#mpl.use('Agg')
import os
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

import numpy as np
from obspy.core import read, Stream, UTCDateTime

#=================================================================================
#                       files and parameters
#=================================================================================
mseed_dir= '/media/tgoebel/Data/seis/BlueMountain/mseed'
plot_dir = f"{os.environ['HOME']}/code/python/Github/EqLocGit/plots"
l_Sta    = [ 'BM%02d'%(i) for i in range(1,11)]

print( l_Sta)
new_sampling = 100 # downsample
# bandpass filter
f_min, f_max = 1, 50
s_com     = 'all' #'Z'  # N, Z, E
# set dCha == None or False to use s_com
dCha      = { 'BM01' : 'Z', 'BM02' : 'N',  'BM03' : 'N', 'BM04' : 'Z',
              'BM06' : 'Z', 'BM07' : 'Z', 'BM09' : 'Z', 'BM10' : 'E'}
save_fig  = False
#tmin     = UTCDateTime('2019-05-02')
tmin     = datetime( 2017,  9, 11)
tmax     = datetime( 2017,  9, 25)

curr_t    = tmin

n_noReads = 0
#===============================2=================================================
#                    load mseed files, stack stream objects
#=================================================================================
while (tmax-curr_t).days > 0:
    t_str = curr_t.strftime( '%Y%m%d')
    print( t_str, 'dt in days', (tmax-curr_t).days)
    st_all = Stream()
    l_Sta_withWf = []
    for sta in l_Sta:
        print(  '--------------------- station', sta, '----------------------------')
        if os.path.isdir( f"{mseed_dir}/{sta}"):
            #write to mseed
            mseed_file = '%s/%s/%s_%s.mseed'%(mseed_dir, sta, sta, t_str)
            if os.path.isfile( mseed_file):
                SEIS = read(  mseed_file)
                # select components
                if dCha is not None and dCha != False:
                    SEIS    = SEIS.select( component= dCha[sta])
                elif s_com is not None and s_com != 'all':
                    SEIS    = SEIS.select( component= s_com)
                print( 'wf file: ', mseed_file.split('/')[-1],'nCha', len(SEIS), 'df', SEIS[0].stats.sampling_rate)
                ###A### down -sample
                if new_sampling is not None and new_sampling != False:
                    SEIS.resample( new_sampling)
                    #for i_s in range(3):
                    #    SEIS[0].decimate( factor=1, no_filter=True)
                    #tr_new.decimate(factor=2.5, strict_length=False)
                    print( 'finished down -sample', SEIS[0].stats.sampling_rate)
                ###B### detrend and filtering, bandpass
                SEIS.detrend( type = 'linear')# linear or simple or demean
                if f_min is not None:
                    print( f"bandpass: {f_min}-{f_max}")
                    SEIS.filter( 'bandpass', freqmin=f_min, freqmax=f_max, corners = 2)
                #print( SEIS[0].stats)

                st_all += SEIS
                #SEIS[0].plot()
                [l_Sta_withWf.append( sta) for i in range( len(SEIS))]
                #l_Sta_withWf.append( sta)
    print( 'tot no. of cha: ', len( st_all))
    #===============================3=================================================
    #                   plot every component
    #=================================================================================
    #a_tr1 = SEIS[0]
    n_axis = len( st_all)
    for i_ax in range( n_axis):
        dt    = st_all[i_ax].stats.delta
        npts  = st_all[i_ax].stats.npts
        sta   = l_Sta_withWf[i_ax]#st_all[i_ax].stats.station
        cha   = st_all[i_ax].stats.channel
        #-------check if wf are shifted because of missing data
        t_start   = st_all[i_ax].stats.starttime
        shift_sec = t_start - UTCDateTime( curr_t.strftime( '%Y-%m-%d') )
        print('---plot-----', sta, cha, npts, dt, '-',
              st_all[i_ax].stats.starttime, st_all[i_ax].stats.endtime,'shift (hr)', shift_sec/3600)
        a_t    = np.arange( 0, dt*npts, dt) + shift_sec #/3600.

        if i_ax == 0:
            plt.figure(1, figsize=(10,14))
            ax = plt.subplot( n_axis, 1, i_ax+1)#, sharex = ax)
        else:
            ax = plt.subplot( n_axis, 1, i_ax+1, sharex = ax)
        ax.plot( a_t, st_all[i_ax][0:npts], 'k-', label = f"{sta}-{cha}")
        ax.legend( loc = 'upper right', fontsize = 'xx-small', frameon = False)
        ax.set_xlabel( f"seconds after {curr_t.strftime('%Y-%m-%d')}, bp={f_min}-{f_max}")
    # n_axis = int( a_t_hr[-1])+1
    # # plot in hourly axis increments
    # #plt.figure(1, figsize=(6,14))
    # for i_ax in range( n_axis):
    #     #ax = plt.subplot( n_axis, 1, i_ax+1)#, sharex = ax)
    #     plt.figure()
    #     ax = plt.subplot( 111)
    #     i1,i2 = int(i_ax*(3600./dt)), int((i_ax+1)*3600./dt)
    #     ax.plot( a_t_hr[i1:i2], a_tr1.data[i1:i2])
    #     ax.set_xlabel( tmin)
    if save_fig == True:
        print( f"save: {plot_dir}/{curr_t.strftime('%Y%m%d')}_{s_com}.png")
        plt.savefig(f"{plot_dir}/{curr_t.strftime('%Y%m%d')}_{s_com}.png")
        plt.clf()
    else:
        plt.show()
    curr_t = curr_t + timedelta(days=1)



