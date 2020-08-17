#!/usr/bin/python3.7
#=================================================================================
#                       modules
#=================================================================================
#import sys
#import argparse
#import logging
#mpl.use('Agg')
import glob, os
#-----------------------------my modules------------------------------------------
from EqLocGit.CD11.src.cd11_reader import CD11

#=================================================================================
#                       files and parameters
#=================================================================================
#CD11_dir = '20190529'
CD11_dir = '/media/tgoebel/Data/seis/BlueMountain/CD11/BM01/20190529' #20180516'
plot_dir = '%s/projects/indSeism/BlueMountain/plots/seis'%( os.path.expanduser('~'))
# bandpass filter
f_min, f_max = .5, 50
#filename = "20190529/20190529_000000_C_BM01.cd11"
l_files = sorted( glob.glob( '%s/*.cd11'%(CD11_dir)))
print( 'total number of files in %s: %i'%( CD11_dir, len( l_files)))
#
n_noReads = 0
for curr_file in l_files:

    # read CD1.1 file format
    cd1_1 = CD11( order="<", check_order=True)#order default="<"
    st = cd1_1.read( curr_file)
    # if OK then plot
    if len(st) > 0:
        print( st[0].stats)
        ###A### spectrogram
        st[0].spectrogram( log = True, title= '%s.%s - %s'%(st[0].stats.station, st[0].stats.channel,st[0].stats.starttime),
                                   outfile = '%s/%s'%( plot_dir, curr_file.split('/')[1].replace('.cd11', 'spec.png')))

        ###B### detrend and filtering, bandpass
        st.detrend( type = 'linear')# linear or simple or demean
        st.filter( 'bandpass', freqmin=f_min, freqmax=f_max, corners = 2)
        ##plot all three traces
        #st.plot( outfile = '%s/%s'%( plot_dir, curr_file.split('/')[1].replace('cd11', 'png')))
        st.plot( title = 'BM1 - %s, fr=%.2f-%i'%(st[0].stats.starttime, f_min, f_max))

    else:
        print( ' %s could not read waveform'%( curr_file))
        n_noReads += 1
print( 'total number of not readable waveforms', n_noReads)

"""
def createParser(usage=False):
    parser = argparse.ArgumentParser(prog='dc11_reader')
    parser.add_argument("-v", "--verbosity", nargs='?', help="Verbosity", type=int, default=0)
    parser.add_argument("-o", "--order", default="<", help="File endian")
    parser.add_argument("-f", "--file", help="input CD11 file name")
    if usage:
        parser.print_help()
    return parser.parse_args()


if __name__ == '__main__':
    usage = len(sys.argv) < 2
    par_args = createParser(usage)
    if usage:
        sys.exit()
    debug = par_args.verbosity
    if debug == 4:
        logging.basicConfig(format='%(levelname)s:%(asctime)s %(message)s', level=logging.DEBUG)
    elif debug == 0:
        logging.basicConfig(format='%(levelname)s:%(asctime)s %(message)s', level=logging.ERROR)
    else:
        logging.basicConfig(format='%(levelname)s:%(asctime)s %(message)s', level=logging.INFO)
    # filename = "/SD/CD1.1/AMAZ.20190918.0000.11bin"
    filename = "/20190529/20190529_000000_C_BM01.cd11"
    #
    # read CD1.1 file format
    cd1_1 = CD11(order=par_args.order, check_order=True)
    st = cd1_1.read(par_args.file)
    # if OK then plot
    if len(st) > 0:
        st.plot()
"""