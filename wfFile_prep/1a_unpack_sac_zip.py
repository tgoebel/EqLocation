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
from datetime import datetime
from datetime import timedelta

import numpy as np

from obspy.core import read
from obspy import Trace, Stream
#=================================================================================================
#                           params
#=================================================================================================
s_end_date = '20170816'
#s_end_date = '20171114' #include first shut-down
dir_in  = f"/media/tgoebel/My Passport/BMSAC_{s_end_date}"
# for os.system and command line
dir_in_bash = f"/media/tgoebel/My\ Passport/BMSAC_{s_end_date}"
# extract sac files here
dir_out      = '/media/tgoebel/Data/seis/BlueMountain/mseed'
dir_out_bash = dir_out
# dir_out     = r"/media/tgoebel/Toshiba 1TB/mseed"
# dir_out_bash= r"/media/tgoebel/Toshiba\ 1TB/mseed"
l_sta   = [ 'BM%02d'%(i) for i in range(1, 11)]
l_sta   = [ 'BM%02d'%(i) for i in range(1, 10)]

test_plot = False
clean_up  = True

l_zip_files = glob.glob( f"{dir_in}/*.zip")
print( l_zip_files)

for sta in l_sta:
    zd = dir_in.split('_')[1]
    new_date = f"{zd[4:6]}{zd[6:8]}{zd[2:4]}"
    curr_zip_file = f"{sta}_{new_date}.zip"
    print( dir_in, curr_zip_file)

    # ==============================1==================================================================
    #                         extract all zip files in current dir
    # =================================================================================================
    #if os.path.isfile( f"{dir_in}/{curr_zip_file}"):
    #   os.system(f"unzip {dir_in_bash}/{curr_zip_file} -d {dir_out_bash}")

    #==============================2==================================================================
    #                     complete files list, start and end dates
    #=================================================================================================
    if os.path.isdir( f"{dir_out}/{sta}/cont/sac"):# check if sac dir exist
        cha = 'EHZ'
        l_wf_files_z = sorted( glob.glob( f"{dir_out}/{sta}/cont/sac/*_{cha}.sac"))
        if len( l_wf_files_z) == 0:
            cha = 'DPZ'
            l_wf_files_z = sorted( glob.glob( f"{dir_out}/{sta}/cont/sac/*_{cha}.sac"))
        if len( l_wf_files_z) == 0:
            error_str = f"neither EHZ nor DPZ sac files in {dir_out}/{sta}/cont/sac"
            raise ValueError( error_str)
        print( 'total no. of sac files: ', 3*len(l_wf_files_z), l_wf_files_z[0:2],l_wf_files_z[-2::])
        # get start, end date and number of days
        dt_start = l_wf_files_z[0].split('/')[-1].split( '_')[0]
        dt_start = datetime( int(dt_start[0:4]), int(dt_start[4:6]), int(dt_start[6:8]))
        dt_end   = l_wf_files_z[-1].split('/')[-1].split( '_')[0]
        dt_end   = datetime( int(dt_end[0:4]), int(dt_end[4:6]), int(dt_end[6:8]))
        no_days  = (dt_end-dt_start).days
        #print( dt_start.strftime('%Y%m%d'), dt_start + timedelta(days=1), dt_end)
        print( 'no. days in current zip', no_days)
        #==============================3==================================================================
        #          file list for each day, three components --> create obspy trace objects
        #=================================================================================================
        curr_dt = dt_start
        for i in range( no_days):
            l_wf_file_curr_day = []
            t_str = '%04d%02d%02d' % (curr_dt.year, curr_dt.month, curr_dt.day)
            dir_file_out = '%s/%s/%s_%s.mseed'%( dir_out, sta, sta, t_str)
            if os.path.isfile( dir_file_out):
                print( dir_file_out, 'already exists')
                curr_dt = curr_dt + timedelta( days = 1)
            else:
                #print( f"{dir_out}/{sta}/{curr_dt.strftime('%Y%m%d')}_*_{sta}_{sta}_EHZ.sac")
                #l_wf_file_curr_day = sorted(glob.glob(f"{dir_out}/{sta}/cont/sac/{curr_dt.strftime('%Y%m%d')}_*_{sta}_{sta}_{cha}.sac"))
                l_wf_file_curr_day = sorted(glob.glob(f"{dir_out}/{sta}/cont/sac/{curr_dt.strftime('%Y%m%d')}_*_{cha}.sac"))
                print( 'curr date', curr_dt, 'no. of files', len(l_wf_file_curr_day))
                i_fi = 0
                for wf_file_z in l_wf_file_curr_day:
                    # read all three components
                    ST_Z = read( wf_file_z)

                    if cha == 'EHZ':
                        ST_E = read( wf_file_z.replace( 'EHZ', 'EHE'))
                        ST_N = read( wf_file_z.replace( 'EHZ', 'EHN'))
                    else:
                        ST_E = read( wf_file_z.replace( 'DPZ', 'DPE'))
                        ST_N = read( wf_file_z.replace( 'DPZ', 'DPN'))
                        ST_Z[0].stats.channel = 'EHZ'
                        ST_N[0].stats.channel = 'EHN'
                        ST_E[0].stats.channel = 'EHE'
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
                curr_dt = curr_dt + timedelta( days=1)
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
                # check if sta dir existis in mseed
                if os.path.isdir( f"{dir_out}/{sta}"):
                    pass
                else:
                    os.path.mkdir(f"{dir_out}/{sta}")
                #write to mseed
                print( 'save mseed: ', dir_file_out, 'from sac: ', wf_file_z)

                SEIS.write( dir_file_out, format='MSEED')

        #===============================5=================================================
        #                    clean up
        #=================================================================================
        # remove all extracted folders before going on to the next .zip file
        if clean_up == True:
            for curr_dir in ['cal', 'cont', 'event', 'log', 'soh', 'status', 'window']:
                os.system( f"rm -r {dir_out}/{sta}/{curr_dir}")


