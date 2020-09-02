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

test_plot = False
clean_up  = True

#l_zip_files = glob.glob( f"{dir_in}/*.zip")
l_zip_files = ['/media/tgoebel/My\ Passport/BMSAC_20170816/BM06_BM07_081617.zip',
               '/media/tgoebel/My\ Passport/BMSAC_20170816/BM03_BM04_081617.zip',
               '/media/tgoebel/My\ Passport/BMSAC_20170816/BM01_BM02_081617.zip']
for curr_zip_file in l_zip_files:
    print( curr_zip_file)
    #if os.path.isfile( curr_zip_file):
    os.system(f"unzip {curr_zip_file} -d {dir_out_bash}")

