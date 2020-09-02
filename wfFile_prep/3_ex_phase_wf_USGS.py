#!/opt/anaconda3/envs/obspy/bin/python (python3.7)
"""
    - event and phase data download example
    using the USGS fdsn client

    @author thgoebel U of Memphis
"""
import os
import numpy as np
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client
import scipy.io
#-------------------

#=============================0==========================================
#                  parameters, files and directories
#========================================================================
dPar   = {   'client'            : "USGS",
             'tmin' : UTCDateTime( 2012, 10, 13),
             'tmax' : UTCDateTime( 2015, 1, 1),
             'xmin' : -119.687, 'xmax' : -119.680,
             'ymin' : 39.667, 'ymax' : 39.672,
             'minMag' : .5,
           }

#=============================1==========================================
#                      download catalog
#========================================================================

client = Client( dPar['client'])
cat = client.get_events(starttime=dPar['tmin'], endtime=dPar['tmax'],
                        minlongitude=dPar['xmin'], maxlongitude=dPar['xmax'],
                        minlatitude=dPar['ymin'], maxlatitude=dPar['ymax'],
                        minmagnitude=dPar['minMag'], #catalog="ISC",
                        includearrivals = True)
print( len( cat))
