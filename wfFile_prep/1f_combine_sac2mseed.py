#python3.7
"""
    Note:
        - there seems to be a buck in obspy.trace,
        tr.__add__          L799 to 850:

        # merge traces depending on NumPy array type
        if True in [isinstance(_i, np.ma.masked_array) for _i in data]:
            data = np.ma.concatenate(data)
        else:
            data = np.concatenate(data)
        #!> enforce original array type in both instances <!

        - un-indent the following line (L841)
        data = np.require(data, dtype=lt.data.dtype)

"""
import matplotlib.pyplot as plt
import os, glob
from datetime import datetime
from datetime import timedelta

import numpy as np

from obspy.core import read
from obspy import Trace, Stream
#=================================================================================================
#                           params
#=================================================================================================
s_date = '20171115'
s_sta  = 'BM02'
dir_in  = f"{os.environ['HOME']}/Desktop/BM/{s_sta}_{s_date}"
# extract sac files here
dir_out      = dir_in

test_plot = True

l_files = sorted( glob.glob( f"{dir_in}/{s_date}_*_EHZ.sac"))
print( 'no. sac files', len(l_files)*3, l_files[0:2], l_files[-2::])
s_sta = l_files[0].split('_')[-2]
print( s_date, s_sta)
#==============================1==================================================================
#                     complete files list, start and end dates
#=================================================================================================
i_fi = 0
for wf_file_z in l_files:
    # read all three components
    ST_Z = read( wf_file_z)
    ST_E = read( wf_file_z.replace( 'EHZ', 'EHE'))
    ST_N = read( wf_file_z.replace( 'EHZ', 'EHN'))

    # add hourly or half hour data to each trace
    if i_fi == 0:
        t0 = ST_Z[0].stats.starttime
        a_trZ = ST_Z[0].copy()
        a_trE = ST_E[0].copy()
        a_trN = ST_N[0].copy()
    else: # add to existing trace objects
        a_trZ += ST_Z[0]
        a_trN += ST_N[0]
        a_trE += ST_E[0]
    i_fi += 1

#---------check for data gaps-----------------
for curr_tr in [a_trZ, a_trN, a_trE]:
    if isinstance(curr_tr.data, np.ma.masked_array):
        curr_tr.data.filled( fill_value = 0)
        curr_tr.data = curr_tr.data.compressed()
#===============================4=================================================
#                     test plot: 24hr of EHZ records
#=================================================================================
print( f"check time stamp: {t0}, {a_trZ.stats.starttime}, endtime: {a_trZ.stats.endtime}")
if test_plot == True:
    n  = a_trZ.stats.npts
    dt = a_trZ.stats.delta
    a_t_hr = np.arange( 0, dt*n, dt)#/3600.
    fig, ax =  plt.subplots( 3,1, sharex=True)
    ax[0].plot( a_t_hr , a_trZ.data[0:n], label = 'Z')
    ax[0].legend( loc = 'upper right')
    ax[1].plot( a_t_hr , a_trN.data[0:n], label = 'N')
    ax[1].legend( loc = 'upper right')
    ax[2].plot( a_t_hr , a_trE.data[0:n], label = 'E')
    ax[2].legend( loc = 'upper right')
    # n_axis = int( a_t_hr[-1])+1
    # # plot in hourly axis increments
    # plt.figure(1, figsize=(6,14))
    # for i_ax in range( n_axis):
    #     ax = plt.subplot( n_axis, 1, i_ax+1)#, sharex = ax)
    #     i1,i2 = int(i_ax*(3600./dt)), int((i_ax+1)*3600./dt)
    #     ax.plot( a_t_hr[i1:i2], a_trZ.data[i1:i2])
    ax[2].set_xlabel( f"seconds after {t0}")
    plt.show()

#===============================5=================================================
#              merge trace object and write stream object to mseed
#=================================================================================
# merge all trace objects
SEIS = Stream( traces=[a_trZ, a_trN, a_trE])
dir_file_out = f"{dir_out}/{s_sta}_{s_date}.mseed"
#write to mseed
print( 'save mseed: ', dir_file_out, 'from sac: ', wf_file_z)
SEIS.write( dir_file_out, format='MSEED')




