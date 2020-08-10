'''
Created on Apr 4, 2013



@author: tgoebel
'''
import os
import numpy as np
#--------------------------------------
from obspy.core.event import Catalog, Event, Pick, WaveformStreamID
from obspy.core import UTCDateTime
#----------------------------------
#import dataIO.dictionaryIO as dicIO
#========================================================================
#                       basic processing
#========================================================================
def get_tr_id( tr):
    """

    :return: {net}.{sta}.{loc}.{cha}
    """
    net, sta, loc, cha = tr.stats.network, tr.stats.station, tr.stats.location, tr.stats.channel
    return f"{net}.{sta}.{loc}.{cha}"

# def prepPickTable( dPick, a_sta, **kwargs):
#     """
#         - get pick table ready so it contains
#           only station from the current inventory
#           given by a_sta
#
#     :param dPick:   - py dic
#                      tags: 'pol', 'sta', 'net', 'cha', 'UTC' etc.
#     :param l_sta:   - python array with strings
#                     - current station inventory
#     :return:  dPick - python dic.
#                     - only current stations in l_sta
#     """
#     dPick['UTC'] = np.array([UTCDateTime(i) for i in dPick['UTC']])
#     # remove spaces for sta
#     for tag in ['sta', 'pol']:
#         dPick[tag] = np.array([i.strip(' ') for i in dPick[tag]])
#     # filter station inventory
#     #tmp_sta, a_id1, a_id2 = np.intersect1d( dPick['sta'], a_sta, return_indices=True)
#     if a_sta is not None:
#         sel = np.isin(dPick['sta'], a_sta, assume_unique=False, invert=False)
#         dPick = dicIO.selDicAll( dPick, sel, testForArray=False)
#     if 'onlyP' in kwargs.keys() and kwargs['onlyP'] == True:
#         dPick = dicIO.selectData( dPick, 'P', 'pha', testForArray= False)
#     return dPick

#========================================================================
#                       phase arrivals
#========================================================================
def aicPicker( vX, **kwargs):
    """ Autoregressive AIC picker.
    see:    H. Zhang, C. Thurber, C. Rowe, Automatic P-wave arrival
            detection BSSA (2003)
            Jochen H. Kurz, Christian U. Grosse, Hans-Wolf Reinhardt
            Strategies for reliable automatic onset time picking of
            acoustic emissions and of ultrasound signals in concrete (2005)
    input: vX = time series

    returns AIC function corresponding to the input sequence X
    AIC(k) = k*log{var(x[1,k]} + (N - k - 1)*log{var(x[k+1,N])}
    """
    
    
    ###compute means and variances from 0 - k
    vX2 = vX**2
    vX_mean  = vX.cumsum()/np.arange(1, vX.shape[0]+1,1)
    vX2_mean = vX2.cumsum()/np.arange(1, vX2.shape[0]+1,1)
    
    # compute variances for each sample from 0 - k
    vVar_1k = vX2_mean - vX_mean**2
    
    ### compute variance again from k - N (end of signal)
    vX  = vX[::-1]
    vX2 = vX**2
    vX_mean  = vX.cumsum()/np.arange(1, vX.shape[0]+1,1)
    vX2_mean = vX2.cumsum()/np.arange(1, vX2.shape[0]+1,1)
    vVar_kN = (vX2_mean - vX_mean**2)[::-1]
    
    # compute N - k -1
    vk   = np.arange(1,vX.shape[0]+1)
    vNk1 = vX.shape[0] - vk - 1
    vk[1]    = 1 # correction for the log10 'domain' error.
    vNk1[-1] = 1 # correction for the log10 'domain' error.
    
    # AIC output vector (may be a bit inaccurate!)
    vAIC = vk*np.log10(vVar_1k) + vNk1*np.log10(vVar_kN)
    
    #Correct for the first and last sample in output sequence
    #(logarithm domain error because of zero variance for 1 sample 
    # input sequence.
    vAIC[1]  = vAIC[2]
    vAIC[0]  = vAIC[1]
    vAIC[-2] = vAIC[-3]
    vAIC[-1] = vAIC[-2]
    return vAIC

#========================================================================
#                        data I/O
#========================================================================
def pick2obspy( lTrig, saveNLLoc = True):
    """
            - create obspy event objects based on info in lTrig
            1) create Event Object with attribute picks
            2) if saveNLLoc save to .obs ASCII using
               event.write( file_out, format = 'NLLOC_OBS")
            ! note that phase files will be saved in cwd,
              use os.chdir() to adjust before calling this fct.
    :param lTrig:     type = list
                      - list of dictionary with individual event picks
    :param saveNLLoc: - boolean, default = True

    :return: catalog - type obspy catalog object
                     - contains only pick info no mags or origin
    """
    n_ev = len( lTrig)
    catalog = Catalog()
    ## event ID and origin time
    for i_ev in range( n_ev):
        oriTime = lTrig[i_ev]['time']
        event   = Event()
        ## pick information
        n_pick  = len( lTrig[i_ev]['stations'])
        for i_pick in range(n_pick):
            picks = Pick()
            picks.time        = oriTime + lTrig[i_ev]['P'][i_pick]
            picks.time_errors = 1e-2
            picks.phase_hint  = lTrig[i_ev]['phase'][i_pick]
            picks.waveform_id = WaveformStreamID(seed_string=lTrig[i_ev]['trace_ids'][i_pick])
            event.picks.append(picks)

        if saveNLLoc == True:
            s_datetime = oriTime.strftime("%Y%m%d.%H%M%S")
            file_out = f"{i_ev+1}_{s_datetime}.obs"
            print(f"save: {os.getcwd()}/{file_out}")
            event.write( file_out, format = 'NLLOC_OBS')
        catalog.append( event)
    return catalog

def pha2obsEv( d_pha, save_NLLoc = False, **kwargs):
    """
            - create obspy event object from info in phase
              dictionary - d_pha

    :param d_pha:     type = dic
                      { 'UTC' : np.array, pick infor
                        'sta' : np.array, str vector with station ID
                        'pha' : np.array, P or S
                        'evID': np.array, preferred ID
                      }
    :param saveNLLoc: - boolean,
                        default = False
                      - phase data is saved in cwd as:
                        [evID].obs

    :return: catalog - type obspy catalog object
                     - contains only pick info no mags or origin
    """
    a_evID = np.unique( d_pha['evID'])
    n_ev = len( a_evID)
    catalog = Catalog()
    ## event ID and origin time
    for i_ev in range( n_ev):
        #oriTime = lTrig[i_ev]['time']
        event   = Event()
        event.resource_id = str( a_evID[i_ev])
        ## event and picks
        sel_ev = a_evID[i_ev] == d_pha['evID']
        n_pick  = sel_ev.sum()
        for i_pick in range(n_pick):
            picks = Pick()
            picks.time        = d_pha['UTC'][sel_ev][i_pick]
            picks.time_errors = 1e-2
            picks.phase_hint  = d_pha['pha'][sel_ev][i_pick]
            picks.waveform_id = WaveformStreamID(station_code = d_pha['sta'][sel_ev][i_pick],
                                                 channel_code = d_pha['cha'][sel_ev][i_pick],
                                                 network_code = d_pha['net'][sel_ev][i_pick],
                                                 location_code= d_pha['loc'][sel_ev][i_pick],
                                                 )
            picks.polarity    = d_pha['pol'][sel_ev][i_pick]
            event.picks.append(picks)

        if save_NLLoc == True:
            if 'catName' in kwargs.keys() and kwargs['catName'] is not None:
                file_out = f"{kwargs['catName']}_{a_evID[i_ev]}.obs"
            else:
                file_out = f"{a_evID[i_ev]}.obs"
            print(f"save: {os.getcwd()}/{file_out}")
            event.write( file_out, format = 'NLLOC_OBS')
        catalog.append( event)
    return catalog




    