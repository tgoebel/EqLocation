"""
    - rsync files on local /media/tgoebel/Data
      and
       /gaia/data/mexico/GFZ/
    rsync -azP --delete --dry-run /media/tgoebel/Data/seis/BlueMountain/mseed/[sta]
                                  thgoebel@ceriml20.ceri.memphis.edu:/gaia/data/mexico/mseed/BM/[sta]

"""
import os, glob
from datetime import datetime, timedelta

#=================================================================================================
#                           specify stations and dates
#=================================================================================================
t0    = datetime( 2017, 9, 4)
nDays = 28 #copy nDays of data from t0
l_sta = [ 'BM%02d'%(i) for i in range(1, 11)]
file_format = 'mseed'
#=================================================================================================
#                           copy files
#=================================================================================================
for sta in l_sta:
    dir_from  = f"/media/tgoebel/Data/seis/BlueMountain/{file_format}/{sta}"
    dir_to    = f"thgoebel@ceriml25.ceri.memphis.edu:/gaia/data/mexico/mseed/BM/{sta}"
    for i_t in range( nDays):
        curr_t = t0 + timedelta( days = i_t)
        curr_file = f"{sta}_{curr_t.strftime('%Y%m%d')}.{file_format}"
        if os.path.isfile( f"{dir_from}/{curr_file}"):
            #print( f"{dir_from}/{curr_file}")
            #print( f"{dir_to}/{curr_file}")
            #os.system( f"rsync -azP --delete --dry-run {dir_from}/{curr_file} {dir_to}/{curr_file}")
            os.system( f"rsync -azP --delete {dir_from}/{curr_file} {dir_to}/{curr_file}")










