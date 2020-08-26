#python3.7
"""
    - make a list of local file / folders for all available stations
    - right now data is saved on second hard-drive use "data_dir" to adjust
    - dir structure:
        seis/BlueMountain/[fileType]/[sta]/[date]

    output:
        create ASCII file with following columns

        date    BM01        BM02   etc.
        [date]  [dir_size]  [dir_size]

"""
import matplotlib.pyplot as plt
import os
import datetime
import numpy as np


def get_size(start_path = '.'):
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(start_path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            # skip if it is symbolic link
            if not os.path.islink(fp):
                total_size += os.path.getsize(fp)

    return total_size

#=================================================================================================
#                           params
#=================================================================================================
dir_main = f"/media/tgoebel/Data/seis/BlueMountain"
dir_out  = f"{os.environ['HOME']}/projects/indSeism/BlueMountain/data/seis"
fileType = 'CD11' # sac CD11 or mseed
l_sta = [ 'BM%02d'%(i) for i in range(1, 11)]
print( l_sta)
#==============================0==================================================================
#                         create time vector
#=================================================================================================
if fileType == 'CD11':
    t1 = datetime.datetime(2018, 5, 1)
    t2 = datetime.datetime(2020, 8, 1)
else:
    t1 = datetime.datetime(2016, 9, 1)
    t2 = datetime.datetime(2016, 11, 1)
file_out = f"BM_file_inventory_{t1.strftime('%Y_%m')}_{t2.strftime('%Y_%m')}_{fileType}.txt"

dt_day = 1
curr_t = t1
delta_sec = (t2-curr_t).total_seconds()

a_date = np.array([], dtype = int)
while delta_sec > 0:
    str_date = curr_t.strftime('%Y%m%d')
    #print( f"------{str_date}---------")
    a_date = np.append( a_date, int( str_date))
    curr_t += datetime.timedelta(days=dt_day)
    delta_sec = (t2-curr_t).total_seconds()
print( 'no. days to check for files: ', len(a_date), a_date[0:5],a_date[-5::])
#==============================1==================================================================
#                           get folders and sizes
#=================================================================================================
mData = np.ones( (len(a_date), len(l_sta)+1), dtype=int)*-1
mData[:,0] = a_date
i = 0
for it in range( len( a_date)):
    for i_sta in range( len(l_sta)):
        sta = l_sta[i_sta]
        ## account for different directory structure for CD11, mseed and sac
        if fileType == 'CD11':
            curr_dir = f"{dir_main}/{fileType}/{sta}/{a_date[it]}"
            if os.path.isdir( curr_dir):
                #dir_size = sum(os.path.getsize(f) for f in os.listdir(curr_dir) if os.path.isfile(f))
                #print( dir_size)
                dir_size = get_size( curr_dir)
                print(curr_dir, dir_size)
                mData[it,i_sta+1] = dir_size
        elif fileType == 'mseed':
            curr_file = f"{dir_main}/{fileType}/{sta}/{sta}_{a_date[it]}.{fileType}"
            print( curr_file)
            if os.path.isfile( curr_file):
                mData[it,i_sta+1] = os.path.getsize( curr_file)
        else: # sac
            print( 'only CD11 and mseed implemented')
            raise ValueError

        #os.path.getsize( curr_dir)
        i += 1
#==============================2==================================================================
#                           write to ASCII
#=================================================================================================
np.savetxt( f"{dir_out}/{file_out}", mData, fmt = '%10i%12i%12i%12i%12i%12i%12i%12i%12i%12i%12i',
            header = '%8s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s'%('date','01','02','03','04','05','06','07','08','09','10',))
print( 'tot. no. of station days',  i)


#==============================3==================================================================
#                           create plot??
#=================================================================================================
plt.figure(1)
ax = plt.axes([.15, .12, .8, .85])
plot = ax.pcolor( mData[:,1::]*1e-6)
plt.colorbar( plot, orientation = 'vertical', label = 'daily record size (MB)')
# shift xtick labels by .5
ax.set_xticks(      np.arange( .5, 10.5))
ax.set_xticklabels( np.arange(  1, 11))
ax.set_xlabel( 'Station ID (BM)')
i_y = np.arange( 0, len(a_date),50)
ax.set_yticks( i_y)
ax.set_yticklabels( a_date[i_y])
plt.savefig(  f"{dir_out}/{file_out}".replace('txt', 'png'))
#ax.set_xticklabels( )
plt.show()









