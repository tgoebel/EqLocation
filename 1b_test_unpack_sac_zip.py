#python3.7
"""
    - unzip all sac files, combine three components to day-long mseed file
    - data structure
    input:
        [main_dir]/sac/[sta]/cont/sac/[datetime]_[sta].sac

    output:
        -daily mseed files:
        [main_dir]/mseed/[sta]_[date].mseed


    Note:
        - there seems to be a buck in obspy.trace, L799 to 804:
        # merge traces depending on NumPy array type
        if True in [isinstance(_i, np.ma.masked_array) for _i in data]:
            data = np.ma.concatenate(data)
        else:
            data = np.concatenate(data)
        #!> enforce original array type in both instances <!
        data = np.require(data, dtype=lt.data.dtype)

"""
import matplotlib.pyplot as plt
import os, glob
#import datetime
import numpy as np

from obspy.core import UTCDateTime, read
from obspy import Trace, Stream
#=================================================================================================
#                           params
#=================================================================================================
#dir_main = f"/media/tgoebel/Data/seis/BlueMountain"
dir_main  = f"{os.environ['HOME']}/projects/indSeism/BlueMountain/data/seis"

#fileType  = 'sac' # sac CD11 or mseed
sta       = 'BM07'
test_file = ["BM07122117_07_oneday.zip"]
test_plot = True
out_format= 'mseed'
clean_up  = True


os.chdir( dir_main)
if os.path.isdir( out_format):# check is mseed dir already exists
    pass
else:
    os.mkdir( out_format)

for curr_zip_file in test_file:
    # ==============================1==================================================================
    #                         extract all zip files in current dir
    # =================================================================================================
    #os.system(f"unzip sac/{curr_zip_file} -d {dir_main}/mseed")

    #==============================2==================================================================
    #                     load hourly sac files combine to daily
    #=================================================================================================
    l_wf_files_z = sorted( glob.glob( f"{out_format}/{sta}/cont/sac/*_EHZ.sac"))
    i_fi = 0
    for wf_file_z in l_wf_files_z:#[18:19]:
        # read all three components
        ST_Z = read( wf_file_z)
        ST_E = read( wf_file_z.replace( 'EHZ', 'EHE'))
        ST_N = read( wf_file_z.replace( 'EHZ', 'EHN'))
        print( ST_Z[0].stats.starttime,ST_Z[0].stats.endtime)
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
    #===============================3=================================================
    #                     test plot: 24hr of EHZ records
    #=================================================================================
    print( f"check time stamp: {t0}, {a_trZ.stats.starttime}, endtime: {a_trZ.stats.endtime}")
    if test_plot == True:
        dt = a_trZ.stats.delta
        a_t_hr = np.arange( 0, dt*a_trZ.stats.npts, dt)/3600.
        fig, ax =  plt.subplots( 3,1, sharex=True)
        ax[0].plot( a_t_hr , a_trZ.data, label = 'Z')
        ax[0].legend( loc = 'upper right')
        ax[1].plot( a_t_hr , a_trN.data, label = 'N')
        ax[1].legend( loc = 'upper right')
        ax[2].plot( a_t_hr , a_trE.data, label = 'E')
        ax[2].legend( loc = 'upper right')
        # n_axis = int( a_t_hr[-1])+1
        # # plot in hourly axis increments
        # plt.figure(1, figsize=(6,14))
        # for i_ax in range( n_axis):
        #     ax = plt.subplot( n_axis, 1, i_ax+1)#, sharex = ax)
        #     i1,i2 = int(i_ax*(3600./dt)), int((i_ax+1)*3600./dt)
        #     ax.plot( a_t_hr[i1:i2], a_trZ.data[i1:i2])
        ax[2].set_xlabel( t0)
        plt.show()

    #===============================4=================================================
    #                    merge trace object and write stream object to mseed
    #=================================================================================
    # merge all trace objects
    SEIS = Stream( traces=[a_trZ, a_trN, a_trE])
    #write to mseed
    t_str = '%04d%02d%02d' % (t0.year, t0.month, t0.day)
    mseed_file = '%s_%s.mseed'%( sta, t_str)
    # check if sta dir existis in mseed
    if os.path.isdir( f"mseed/{sta}"):
        pass
    else:
        os.path.mkdir(f"mseed/{sta}")
    print( 'save mseed: ', '%s/mseed/%s/%s'%( dir_main, sta, mseed_file))

    SEIS.write( 'mseed/%s/%s'%( sta, mseed_file), format='MSEED')

    #===============================5=================================================
    #                    clean up
    #=================================================================================
    # remove all extracted folders before going on to the next .zip file
    if clean_up == True:
        for curr_dir in ['cal', 'cont', 'event', 'log', 'soh', 'status', 'window']:
            os.system( f"rm -r mseed/{sta}/{curr_dir}")





