#python3.7
"""
    - only unpach header files to keep track of station locations

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
s_end_date = '20161017'
#s_end_date = '20171114' #include first shut-down
dir_in  = f"/media/tgoebel/My Passport/BMSAC_{s_end_date}"
# for os.system and command line
dir_in_bash = f"/media/tgoebel/My\ Passport/BMSAC_{s_end_date}"
# extract sac files here
dir_out      = '/media/tgoebel/Data/seis/BlueMountain/mseed'
dir_out_bash = dir_out
# dir_out     = r"/media/tgoebel/Toshiba 1TB/mseed"
# dir_out_bash= r"/media/tgoebel/Toshiba\ 1TB/mseed"
#l_sta   = [ 'BM%02d'%(i) for i in range(1, 11)]
l_sta   = ['BM02', 'BM03', 'BM04', 'BM05', 'BM08', 'BM09', 'BM10']
l_sta_old = ['BM08', 'BM03']
l_sta_new = ['BM09', 'BM10']



#===============================0==================================================================
#                           rename files in 03 and 08
#=================================================================================================
i_s = 0
for sta_new in l_sta_new:
    sta_old = l_sta_old[i_s]
    l_files = glob.glob( f"{dir_out}/{sta_old}/{sta_new}*.mseed")
    print('no file to rename', len( l_files), l_files[0], l_files[-1])
    for curr_file in l_files:
        print( f"{curr_file} {curr_file.replace(sta_new, sta_old)}")
        os.system( f"mv {curr_file} {curr_file.replace(sta_new, sta_old)}")
    i_s += 1
print( grrg)




#=================================================================================================
#
#=================================================================================================
l_zip_files = glob.glob( f"{dir_in}/*.zip")
print( l_zip_files)

for sta in l_sta:
    zd = dir_in.split('_')[1]
    new_date = f"{zd[4:6]}{zd[6:8]}{zd[2:4]}"
    curr_zip_file = f"{sta}_{new_date}.zip"
    print( dir_in, curr_zip_file)

    # ==============================1==================================================================
    #                         extract header files to /Data
    # =================================================================================================
    head_dir = f"{dir_out}/{sta}/headers"
    if os.path.isdir( f"{head_dir}"):
        pass
    else:
        os.mkdir(f"{head_dir}" )
    if os.path.isfile( f"{dir_in}/{curr_zip_file}"):
       os.system(f"unzip {dir_in_bash}/{curr_zip_file} '{sta}/cont/headers/*' -d {head_dir}")
       os.system( f"mv {head_dir}/{sta}/cont/headers/*.txt {head_dir}")
       #os.system( f"rm -r {head_dir}/{sta}")

    #==============================2==================================================================
    #                     get sta location from header
    #=================================================================================================
    if os.path.isdir( head_dir):# check if sac dir exist
        l_files = sorted( glob.glob( f"{head_dir}/*.txt"))
        for curr_file in l_files:
            file_obj = open( curr_file, 'r')
            l_str =  file_obj.readline().split('')
            file_obj.close()
            print( l_str)





