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
import os, glob
import datetime
import numpy as np
#=================================================================================================
#                           params
#=================================================================================================
dir_main = f"/media/tgoebel/Data/seis/BlueMountain"
#dir_wf   = "/media/tgoebel/Data/seis/BlueMountain/sac/[sta]"
fileType = 'sac' # sac CD11 or mseed
l_sta = [ 'BM%02d'%(i) for i in range(1, 11)]
print( l_sta)

#==============================1==================================================================
#                         extract all zip files in current dir
#=================================================================================================



#==============================2==================================================================
#                     load hourly sac files combine to daily
#=================================================================================================

# load three component files


#==============================3==================================================================
#                         save as mseed
#=================================================================================================


"""
mData = np.ones( (len(a_date), len(l_sta)+1), dtype=int)*-1
mData[:,0] = a_date
i = 0
for it in range( len( a_date)):
    for i_sta in range( len(l_sta)):
        sta = l_sta[i_sta]
        curr_dir = f"{dir_main}/{fileType}/{sta}/{a_date[it]}"

        if os.path.isdir( curr_dir):
            #dir_size = sum(os.path.getsize(f) for f in os.listdir(curr_dir) if os.path.isfile(f))
            #print( dir_size)
            dir_size = get_size( curr_dir)
            print(curr_dir, dir_size)
            mData[it,i_sta+1] = dir_size
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
"""








