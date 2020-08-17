"""
    - compare NLLoc locations for April 21
      with location from SRL publication for the same data

    1) get station inventory and plot stations
    2) load NLLoc locations, plot
    3) load published catalog --> plot
    @author thgoebel U of Memphis
"""
import os
import numpy as np
import matplotlib as mpl
#mpl.use( 'Agg')
import matplotlib.pyplot as plt

##
from obspy.core import UTCDateTime
from obspy.core.event import read_events, Magnitude, FocalMechanism, NodalPlanes
#
from obspy.clients.fdsn import Client


def obsCat2seisCat( obsCat):
    """
        - convert obspy catalog object to simple data table with:
        ['Time', 'Lon', 'Lat', 'Depth', 'MAG', 'UTC', 'errh', 'errz', 'stderr']
        return:
                dData = {} - python dictionary
                           - tags = 'Time', 'Lon', 'Lat', 'Depth', 'MAG', 'UTC', 'errh', 'errz', 'stderr'

    """
    N = len(obsCat.events)
    l_head = ['Time', 'Lon', 'Lat', 'Depth', 'MAG', 'UTC', 'errh', 'errz', 'stderr']
    dData = {}
    for tag in l_head:
        dData[tag] = np.zeros(N)
    dData['N'] = np.zeros(N, dtype=int)
    i_ev = 0
    for ev in obsCat.events:
        # only if catname is in title
        try:
            id = int(ev.creation_info['author'].split(' ')[0].split('_')[1].split('.')[0])
        except:
            if ev.creation_info['author'].split(' ')[0].split('/')[-1] == 'synth.obs':
                id = 0
            else:
                try:
                    id = int(ev.creation_info['author'].split(' ')[0].split('/')[-1].strip('.obs'))
                except:
                    id = i_ev + 1
        dData['N'][i_ev] = id
        dData['Lon'][i_ev] = ev.origins[0].longitude
        dData['Lat'][i_ev] = ev.origins[0].latitude
        dData['Depth'][i_ev] = ev.origins[0].depth * 1e-3
        # need conversion to decimal year function
        #dData['Time'][i_ev] = date_utils.UTC2dec(ev.origins[0].time)
        dData['UTC'][i_ev] = ev.origins[0].time
        # check if mag exists
        if len(ev.magnitudes) > 0:
            dData['MAG'][i_ev] = ev.magnitudes[0].mag
        ## hor and vertical error
        if ev.origins[0].longitude_errors.uncertainty is not None:
            dData['stderr'][i_ev] = ev.origins[0].quality.standard_error
            dData['errh'][i_ev] = 110 * (.5 * (ev.origins[0].longitude_errors.uncertainty +
                                                   ev.origins[0].latitude_errors.uncertainty))
        if ev.origins[0].depth_errors.uncertainty is not None:
            dData['errz'][i_ev] = ev.origins[0].depth_errors.uncertainty * 1e-3
        i_ev += 1
    return dData
#=============================0==========================================
#                  parameters, files and directories
#========================================================================
#data_dir   = f"{os.environ['HOME']}/data/seismograms/MtCarmel/"
data_dir   = f"{os.environ['HOME']}/PROG/NLLoc/MtCarmel"
plot_dir   = f"{os.environ['HOME']}/projects/eqLocation/MtCarmel"
# ! eventually this file should be moved to a save location to avoid over-writing!
nlloc_file   = 'MtCarmel.sum.grid0.loc.hyp'
cat_out_file =  None #'MtCarmel_4_21_2008.mat'

tmin         = UTCDateTime( "2008-04-21T00:00:00")
tmax         = tmin + 24*3600         #UTCDateTime( "2008-04-21T01:00:00")
MS_lon, MS_lat = -87.89, 38.45
t_MS         = UTCDateTime(2008, 4, 18, 9, 36, 59)
rmax         = .5 # in degree
net, sta,cha = 'Z3', 'MC*', 'HHZ'
#-----------plotting parameters----------------
max_h_err  = 5 # in km


#=============================1==========================================
#                  station inventory
#========================================================================
print(  tmin, tmax, f't-MS{t_MS}')

client = Client("IRIS")

inventory = client.get_stations(starttime=tmin, endtime=tmax, network=net,
                                longitude=MS_lon, latitude=MS_lat, maxradius=rmax,
                                )
inv = inventory.select(station=sta)

plt.figure(1)
ax = plt.subplot( 111)

for obj_sta in inv.networks[0].stations:
    curr_sta =  obj_sta.get_contents()['stations'][0].split(' ')[0]
    print( '---------------',curr_sta,'-----------------')
    ax.plot( obj_sta.longitude, obj_sta.latitude, '^', mfc = '.5')
    ax.text(obj_sta.longitude, obj_sta.latitude, curr_sta)


#=============================2==========================================
#                  load and plot NLLoc catalog
#========================================================================
catObs = read_events( f"{data_dir}/{nlloc_file}")
print( 'no. of seismic events in obspy cat', len(catObs))
dCat = obsCat2seisCat( catObs)
print( 'no events in py dic: ', len(dCat['N']))

if cat_out_file is not None:
    os.chdir( data_dir)
    import scipy.io
    scipy.io.savemat( cat_out_file, dCat)
    #cat.write( cat_out_file+'.cmt',      format = 'CMTSOlUTION')
    #cat.write(cat_out_file + '.cnv',     format='CNV')
    #cat.write(cat_out_file + '.nordic',  format='NORDIC')
    #cat.write(cat_out_file + '.sc2ml',   format='SC3ML')
    #cat.write(cat_out_file + '.scardec', format='SCARDEC')
    #cat.write(cat_out_file + '.zmap',    format='ZMAP')

norm  = mpl.colors.Normalize(vmin = 0, vmax = max_h_err)
ax.set_title( f"{ UTCDateTime(dCat['UTC'][0]).strftime('%Y/%m/%d')}, no. ev: {len(catObs)}")
plot1 = ax.scatter( dCat['Lon'], dCat['Lat'], c = dCat['errh'], s = 15, norm = norm, cmap = plt.cm.RdYlGn_r)
cbar = plt.colorbar(plot1, norm=norm)
cbar.set_label( 'h-err (km)')
plt.savefig( f"{data_dir}/{nlloc_file.replace('hyp', 'png')}")
plt.show()





