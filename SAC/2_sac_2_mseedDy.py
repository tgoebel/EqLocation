#!/usr/bin/python3
"""
    -load half of hourly data files in sac for all available
    stations in Blue Mountain between tmin and tmax
    - convert to day-long miniseed files
Before running this code:
    - move sac_reader.py to sub-directory call 'src' with file __init__.py
"""
#=================================================================================
#                       modules
#=================================================================================
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob, os
import numpy as np
from obspy.core import UTCDateTime
from obspy import Trace, Stream
from shutil import copyfile
import zipfile
impprt os # os.system( "unzip [filename] -d [target_dir]" 
# from obspy.core import read  read( "[date_time].sac")
#-----------------------------my modules------------------------------------------
from EqLocGit.SAC.src.sac_reader import SAC
#=================================================================================
#                       files and parameters
#=================================================================================
# #Where is the path to new file?
# input_file='BM07122117_07_oneday.zip'  
# # Unzip file
# with zipfile.ZipFile(input_file, 'r') as zip_ref:
#     zip_ref.extractall('/media/tgoebel/Data/seis/BlueMountain/SAC')
#     
# SAC_dir = '/media/tgoebel/Data/seis/BlueMountain/SAC/BM07/cont/sac'
# mseed_dir= '/media/tgoebel/Data/seis/BlueMountain/mseed'
SAC_dir = '/Users/shima/Desktop/BM/BM07/cont/sac'
mseed_dir= 'Users/shima/Desktop/BM/BM07/mseed'

test_plot= True
# l_Sta    = ['BM01', 'BM02', 'BM03', 'BM04', 'BM06', 'BM07', 'BM09', 'BM10']
l_Sta    = ['BM07']
tmin     = UTCDateTime('2017-12-12')
tmax     = UTCDateTime('2017-12-13')

for sta in l_Sta:
    print(  '--------------------- station', sta, '----------------------------')
    curr_t    = tmin
    n_noReads = 0
    #===============================2=================================================
    #                      load half an hourly SAC data files
    #=================================================================================
    while (tmax-curr_t) > 0:
        t_str = '%04d%02d%02d'%(curr_t.year,curr_t.month,curr_t.day)
        wf_dir = '%s/%s/%s'%( SAC_dir, sta,t_str)
        l_Files = sorted( glob.glob( '%s/*.sac'%( wf_dir)))
        print('-----------current folder:', t_str, 'n wf files: %i'%( len( l_Files)))
        i_fi = 0
        for curr_file in l_Files:
            # read CD1.1 file format TODO replace with unpacked sac files
            sac1 = SAC( order="<", check_order=True)#order default="<"
            st    = sac1.read( curr_file)
            print( st[0].stats.starttime,st[0].stats.endtime)
            if len(st) > 0:
                if i_fi == 0:
                    t0 = st[0].stats.starttime
                    a_tr1 = st[0].copy()
                    a_tr2 = st[1].copy()
                    a_tr3 = st[2].copy()
                else: # add to existing trace objects
                    a_tr1 += st[0]
                    a_tr2 += st[1]
                    a_tr3 += st[2]
                i_fi += 1
            else:
                print( ' %s could not read waveform'%( curr_file))
                n_noReads += 1
        #===============================3=================================================
        #                     test plot: 24hr of first trace records
        #=================================================================================
        if test_plot == True:
            dt = a_tr1.stats.delta
            a_t_hr = np.arange( 0, dt*a_tr1.stats.npts, dt)/3600.
            n_axis = int( a_t_hr[-1])+1
            # plot in hourly axis increments
            plt.figure(1, figsize=(6,14))
            for i_ax in range( n_axis):
                ax = plt.subplot( n_axis, 1, i_ax+1)#, sharex = ax)
                i1,i2 = int(i_ax*(3600./dt)), int((i_ax+1)*3600./dt)
                ax.plot( a_t_hr[i1:i2], a_tr1.data[i1:i2])
            ax.set_xlabel( tmin)
            plt.show()

        #===============================4=================================================
        #                    merge trace object and write stream object to mseed
        #=================================================================================
        # merge all trace objects
        SEIS = Stream( traces=[a_tr1, a_tr2, a_tr3])
        #write to mseed
        mseed_file = '%s_%s.mseed'%( sta, t_str)
        print( 'save mseed: ', '%s/%s/%s'%( mseed_dir, sta, mseed_file))
        SEIS.write( '%s/%s/%s'%( mseed_dir, sta, mseed_file), format='MSEED')
        curr_t += 3600*24
    print( 'total number of not readable waveforms', n_noReads)
