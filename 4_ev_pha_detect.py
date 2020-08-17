#!/opt/anaconda3/envs/obspy/bin/python (python3.7)
"""
    - identify events using a network sta/lta detector
    - get picks for every detection
    - save picks as m_pick = { '[evID]_[sta]_[cha]' : float([p_arrival]) }
    ! note that both order of stream objects and trigger detects
      is determined by trace-ids (i.e. controlled by net.sta)

"""
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use( 'Agg')
import os, glob, time
import numpy as np
from pprint import pprint

t1 =  time.time()

##
from obspy.core import read, Stream, UTCDateTime
from obspy.signal.trigger import coincidence_trigger

#----------------------------------------
import seis.pickerUtils as pickUtils
#=============================0==========================================
#                  parameters, files and directories
#========================================================================
s_catalog  = 'MtCarmel'
wf_dir     = f"{s_catalog}/mseed"
plot_dir   = f"{s_catalog}/plots"
phase_dir  = f"{s_catalog}/obs"
tmin         = UTCDateTime( "2008-04-21T00:00:00")
tmax         = tmin + 24*3600 # only one day of seismic data

plot_detect  = True
show_picks   = True
saveNLLoc    = True
fmin, fmax   = 2, 45
#----event detection params
trig_on, trig_off = 5, 1     # only for sta/lta
sta, lta          = .5, 10   # only for sta/lta
trig_sum          = 3 #
min_no_picks      = 4 #no. of stations or channel
# #----------filter pick algorithm
# threshold_1, threshold_2 = 5, 10

## time window picks, relative to ev. detection
t_pick_min, t_pick_max         = 2,3 # 6, 10

a_all_sta = ['MC1', 'MC2', 'MC4', 'MC5']
l_files = sorted( glob.glob( f'{wf_dir}/*.mseed'))
print( f'import following traces: {len(l_files)} - {l_files}')

#=============================1==========================================
#        load vertical components, merge all traces to one stream
#========================================================================
st_Z = Stream()
l_sta= []
for filename in l_files:
    st_tmp = read( filename)
    # detrend, taper
    st_tmp.detrend( type = 'linear')
    # filter - bandpass
    if fmin != None:
        st_tmp.filter('bandpass', freqmin=fmin, freqmax=fmax)  # optional prefiltering
    # only z components
    st_tmp = st_tmp.select( component='Z')
    st_Z += st_tmp
    l_sta.append( st_tmp[0].stats.station)
# get rid of potential data gaps
st_Z.merge( fill_value=0)
print( st_Z)

#=============================2==========================================
#                  network trigger
#========================================================================
# type = 'classicstalta', "recstalta", 'carlstatrig', 'zdetect
trig = coincidence_trigger("recstalta", trig_on, trig_off, st_Z,
                           trig_sum, sta=sta, lta=lta)
nTrig  = len( trig)
print( f"no. of triggers: {nTrig}")

#---------------get meta data-------------------------------
# get start of trace
t0_UTC = st_Z[0].stats.starttime
dt     = st_Z[0].stats.delta # samling rate

nTra   = len( st_Z)

print( f'sta. {l_sta}, no. of traces: {nTra}, {len(l_sta)}, no. Trig: {nTrig}, dt= {dt}')

#=============================3==========================================
#                  phase arrivals at each station
#========================================================================
n_no_picks    = 0
for i_ev in range( nTrig):
    if len( trig[i_ev]['stations']) >= min_no_picks:
        trig[i_ev]['P']     = np.ones(len( trig[i_ev]['stations']))*np.nan
        trig[i_ev]['phase'] = np.empty(len( trig[i_ev]['stations']), dtype =str)
        ## loop over every station
        i_pick = 0 # pick counter, not all stations have picks
        for i_sta in range( len( trig[i_ev]['stations'])):
            s_sta  = trig[i_ev]['stations'][i_sta]
            # trim vertical stream shorter window for FPick
            tr = st_Z.copy().select(station = s_sta)
            tr = tr[0].trim(trig[i_ev]['time'] - t_pick_min, trig[i_ev]['time'] + t_pick_max)

            a_cFct = pickUtils.aicPicker( tr.data)
            p_pick = np.argmin(a_cFct)*dt #a_t[ a_cFct == a_cFct.min()][0]
            print(f"AIC pick - rel pick: {p_pick}")
            trig[i_ev]['P'][i_pick]     = tr.stats.starttime + p_pick - trig[i_ev]['time']
            trig[i_ev]['phase'][i_pick] = 'P'
            if show_picks == True:
                UTC_AT = trig[i_ev]['time']
                s_datetime = UTC_AT.strftime("%Y%m%d.%H%M%S")
                fig, ax = plt.subplots(2, 1, sharex='all')  # , sharey='all')
                a_t = np.linspace(0, tr.stats.npts * dt, tr.stats.npts)
                ax[0].plot(a_t, tr.data, 'k-')
                ax[0].axvline(p_pick, c='#FFA500', label='AIC')
                a_cFct -= a_cFct.mean()
                ax[1].plot( a_t, a_cFct, 'g-')
                ax[1].axvline( p_pick, c = '#FFA500', label = 'AIC')
                ax[1].legend()
                ax[0].set_xlim( p_pick-.5, p_pick+2)
                plt.savefig( f"{plot_dir}/{s_catalog}_{s_datetime}_{s_sta}.png")
                plt.show()
                plt.clf()
            i_pick += 1
        else: # len( l_trig_on_off) == 0:
            n_no_picks += 1

print( 'station with net detect but missing pick: ', n_no_picks)
#============================4===========================================
#                      test plot - waveforms and picks
#========================================================================

#st.plot( type="relative")
a_t = np.arange( 0, st_Z[0].stats.npts*st_Z[0].stats.delta, st_Z[0].stats.delta)
fig, ax = plt.subplots( nTra, 1, sharex = 'all', sharey = 'all')
fig.set_figwidth(20)
cmap = mpl.cm.get_cmap( 'cool') # different color for each event
to_hr = 1 #1./3600

for s_sta in a_all_sta:
    i_ax = a_all_sta.index( s_sta)
    tr_tmp = st_Z.copy().select( station = s_sta)
    ax[i_ax].plot(a_t * to_hr, tr_tmp[0].data, 'k-', label = s_sta)
    ax[i_ax].set_ylim(-50 * 1e3, 50 * 1e3)
    ax[i_ax].legend(loc = 'upper right')

l_trig_fin = []
for i_ev in range( nTrig):
    if len( trig[i_ev]['stations']) >= min_no_picks:
        l_trig_fin.append( trig[i_ev])
        ## loop over every station
        i_pick = 0 # pick counter, not all stations have picks
        for i_sta in range( len( trig[i_ev]['stations'])):
            s_sta  = trig[i_ev]['stations'][i_sta]
            # check if pick exists for this event and trace
            print( i_ev, trig[i_ev]['time'], i_sta, s_sta, trig[i_ev]['P'])
            f_P_pick = trig[i_ev]['time'] - t0_UTC + trig[i_ev]['P'][i_sta]
            i_ax = a_all_sta.index(s_sta)
            ax[i_ax].axvline( f_P_pick*to_hr, c = cmap( i_ev/len(trig)), lw = 1.5, ls = '--')

ax[0].set_title( f"{t0_UTC} - no.ev= {len(l_trig_fin)}")
plt.savefig(f'{plot_dir}/pick_{tmin}_{tmax}.png')
print( f"no.ev= {len(l_trig_fin)}")
plt.show()
pprint( l_trig_fin)
#============================5===========================================
#                save phase picks for every event
#========================================================================
# save NLLoc phase data
os.chdir( phase_dir)
catalog = pickUtils.pick2obspy( l_trig_fin, saveNLLoc=saveNLLoc)

print( 'run time', time.time()-t1)



