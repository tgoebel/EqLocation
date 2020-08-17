#!/usr/bin/python3

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import io

import sys
import argparse

import matplotlib as mpl
mpl.use('Qt5Agg')

import logging
import struct
import numpy as np
from obspy import UTCDateTime, Trace, Stream
from obspy.core import Stats

RANGEBIT = 0x8000
log = logging.getLogger("CD112MSEED")


def get_signed_val(value, sign_bit, mask):
    if (value & sign_bit) != 0:
        bmask = -sign_bit
        value = value | bmask
    return value


def _unpack(diffs, start, comp, offset, bits, groupCode=0):
    y0 = 0
    y1 = 0
    y2 = 0
    y3 = 0
    if bits == 4:
        y0 = comp[offset] >> 4
        y1 = comp[offset] & 0x0000000F
        y2 = comp[offset+1] >> 4
        y3 = comp[offset+1] & 0x0000000F
        offset += 2
    elif bits == 6:
        y0 = (comp[offset] & 0x000000FC) >> 2
        y1 = (comp[offset] & 0x00000003) << 4 | comp[offset + 1] >> 4
        offset += 1
        y2 = (comp[offset] & 0x0000000F) << 2 | comp[offset + 1] >> 6
        offset += 1
        y3 = comp[offset] & 0x0000003F
        offset += 1
    elif bits == 8:
        y0 = comp[offset]
        offset += 1
        y1 = comp[offset]
        offset += 1
        y2 = comp[offset]
        offset += 1
        y3 = comp[offset]
        offset += 1
    elif bits == 10:
        y0 = comp[offset] << 2 | (comp[offset + 1] & 0x000000C0) >> 6
        offset += 1
        y1 = (comp[offset] & 0x0000003F) << 4 | comp[offset + 1] >> 4
        offset += 1
        y2 = (comp[offset] & 0x0000000F) << 6 | (comp[offset + 1] & 0x000000FC) >> 2
        offset += 1
        y3 = (comp[offset] & 0x00000003) << 8 | comp[offset + 1]
        offset += 2
    elif bits == 12:
        y0 = comp[offset] << 4 | comp[offset + 1] >> 4
        offset += 1
        y1 = (comp[offset] & 0x0000000F) << 8 | comp[offset + 1]
        offset += 2
        y2 = comp[offset] << 4 | comp[offset + 1] >> 4
        offset += 1
        y3 = (comp[offset] & 0x0000000F) << 8 | comp[offset + 1]
        offset += 2
    elif bits == 14:
        y0 = comp[offset] << 6 | (comp[offset + 1] & 0x000000FC) >> 2
        offset += 1
        y1 = (comp[offset] & 0x00000003) << 12 | comp[offset + 1] << 4 | comp[offset + 2] >> 4
        offset += 2
        y2 = (comp[offset] & 0x0000000F) << 10 | comp[offset + 1] << 2 | (comp[offset + 2] & 0x000000C0) >> 6
        offset += 2
        y3 = (comp[offset] & 0x0000003F) << 8 | comp[offset + 1]
        offset += 2
    elif bits == 16:
        y0 = comp[offset] << 8 | comp[offset + 1]
        offset += 2
        y1 = comp[offset] << 8 | comp[offset + 1]
        offset += 2
        y2 = comp[offset] << 8 | comp[offset + 1]
        offset += 2
        y3 = comp[offset] << 8 | comp[offset + 1]
        offset += 2
    elif bits == 18:
        y0 = comp[offset] << 10 | comp[offset + 1] << 2 | (comp[offset + 2] & 0x000000C0) >> 6
        offset += 2
        y1 = (comp[offset] & 0x0000003F) << 12 | comp[offset + 1] << 4 | comp[offset + 2] >> 4
        offset += 2
        y2 = (comp[offset] & 0x0000000F) << 14 | comp[offset + 1] << 6 | (comp[offset + 2] & 0x000000FC) >> 2
        offset += 2
        y3 = (comp[offset] & 0x00000003) << 16 | comp[offset + 1] << 8 | comp[offset + 2]
        offset += 3
    elif bits == 20:
        y0 = comp[offset] << 12 | comp[offset + 1] << 4 | comp[offset + 2] >> 4
        offset += 2
        y1 = comp[offset] << 16 | comp[offset + 1] << 8 | comp[offset + 2]
        offset += 3
        y2 = comp[offset] << 12 | comp[offset + 1] << 4 | comp[offset + 2] >> 4
        offset += 2
        y3 = comp[offset] << 16 | comp[offset + 1] << 8 | comp[offset + 2]
        offset += 3
    elif bits == 24:
        y0 = comp[offset] << 16 | comp[offset + 1] << 8 | comp[offset + 2]
        offset += 3
        y1 = comp[offset] << 16 | comp[offset + 1] << 8 | comp[offset + 2]
        offset += 3
        y2 = comp[offset] << 16 | comp[offset + 1] << 8 | comp[offset + 2]
        offset += 3
        y3 = comp[offset] << 16 | comp[offset + 1] << 8 | comp[offset + 2]
        offset += 3
    elif bits == 28:
        y0 = comp[offset] << 20 | comp[offset + 1] << 12 | comp[offset + 2] << 4 | comp[offset + 3] >> 4
        offset += 3
        y1 = (comp[offset] & 0x0000000F) << 24 | comp[offset + 1] << 16 | comp[offset + 2] << 8 | comp[offset + 3]
        offset += 4
        y2 = comp[offset] << 20 | comp[offset + 1] << 12 | comp[offset + 2] << 4 | comp[offset + 3] >> 4
        offset += 3
        y3 = (comp[offset] & 0x0000000F) << 24 | comp[offset + 1] << 16 | comp[offset + 2] << 8 | comp[offset + 3]
        offset += 4
    elif bits == 32:
        y0 = comp[offset] << 24 | comp[offset+1] << 16 | comp[offset+2] << 8 | comp[offset+3]
        offset += 4
        y1 = comp[offset] << 24 | comp[offset+1] << 16 | comp[offset+2] << 8 | comp[offset+3]
        offset += 4
        y2 = comp[offset] << 24 | comp[offset+1] << 16 | comp[offset+2] << 8 | comp[offset+3]
        offset += 4
        y3 = comp[offset] << 24 | comp[offset+1] << 16 | comp[offset+2] << 8 | comp[offset+3]
        offset += 4
    else:
        log.error("ERROR: decompress")
        return
    sign_bit = int((1 << bits) / 2)
    mask = sign_bit - 1
    y00 = get_signed_val(y0, sign_bit, mask)
    diffs[start + 0] = y00
    y10 = get_signed_val(y1, sign_bit, mask)
    diffs[start + 1] = y10
    y20 = get_signed_val(y2, sign_bit, mask)
    diffs[start + 2] = y20
    y30 = get_signed_val(y3, sign_bit, mask)
    diffs[start + 3] = y30


def uncompress(comp, nsamp, next_sample, endian):
    if comp is None:
        return None
    if nsamp % 20 != 0:
        return None
    if len(comp) == 0 or nsamp == 0:
        return []
    samples = np.zeros(nsamp)
    ngroup = int(nsamp / 20)
    offset = int(nsamp / 10)
    # currentSample, = struct.unpack_from('>I', comp, offset)
    current_sample = comp[offset] << 24 | comp[offset+1] << 16 | comp[offset+2] << 8 | comp[offset+3]
    offset += 4
    for igroup in range(ngroup):
        blockStart = int(5 * igroup)
        blockEnd = blockStart + 5
        # groupCode, = struct.unpack_from('>H', comp, int(2 * igroup))
        ind = int(2 * igroup)
        groupCode = comp[ind] << 8 | comp[ind + 1]
        highRange = groupCode & RANGEBIT != 0
        divisor = 2
        if highRange:
            divisor = 4
        for iblock in range(blockStart, blockEnd):
            shift = 12 - int(3 * (iblock - blockStart))
            codeBits = (groupCode >> shift) & 0x00000007
            packSize = codeBits * divisor + 4
            _unpack(samples, int(iblock * 4), comp, offset, packSize)
            offset += int(packSize / 2)
    for ix in range(1, nsamp):
         samples[ix] = samples[ix] + samples[ix - 1]
    for ix in range(nsamp):
        save = samples[ix]
        samples[ix] = current_sample
        current_sample += save
    if next_sample is not None and len(next_sample) > 0:
        next_sample[0] = current_sample
    return samples


class FrameTrailer(object):
    endian = ">"
    authKeyId = 0
    authLen = 0
    authValue = 0
    crc64 = 0

    def __init__(self, order=">"):
        self.endian = order

    def read(self, dis):
        try:
            self.authKeyId, = struct.unpack(self.endian + 'I', dis.read(4))
            self.authLen, = struct.unpack(self.endian + 'I', dis.read(4))
            if self.authLen > 40:
                log.error("Corrupt FrameTrailer: {} auth size?".format(self.authLen))
                return False
            fmt = "{0}B".format(self.endian)
            self.authValue = [struct.unpack(fmt, dis.read(1))[0] for p in range(self.authLen)]
            self.crc64, = struct.unpack(self.endian + 'Q', dis.read(8))
        except:
            return False
        return True


class TxState(object):
    endian = ">"
    frame_time_millis = 0
    duration_millis = 0
    num_chan = 0

    def __init__(self, frameTimeMillis, durationMillis, numChan, order=">"):
        self.endian = order
        self.frame_time_millis = frameTimeMillis
        self.duration_millis = durationMillis
        self.num_chan = numChan


class ChannelSubFrame(object):
    endian = ">"
    header = None
    trailer = None
    samples = None
    next_sample = 0
    num_valid_samples = 0
    compressedData = None
    padding = []
    channel_length = 0
    authentication_offset = 0
    authentication = 0
    transformation = 0
    sensor_type = 0
    option_flag = 0
    site_name = None
    channel_name = None
    location_name = None
    data_format = None
    calib_factor = 0
    calib_period = 0
    time_stamp = None
    subframe_time_length = 0
    npts = 0
    channel_status_size = 0
    channel_status_data = None
    data_size = 0
    channel_data = None
    subframe_count = 0
    auth_key_ident = 0
    auth_size = 0
    auth_value = None
    next_value = []
    data = None
    trace = None
    duration = 1

    def __init__(self, order=">"):
        self.endian = order

    def read(self, dis):
        self.channel_length, = struct.unpack(self.endian + 'I', dis.read(4))
        self.authentication_offset, = struct.unpack(self.endian + 'I', dis.read(4))
        self.authentication, = struct.unpack(self.endian + 'B', dis.read(1))
        self.transformation, = struct.unpack(self.endian + 'B', dis.read(1))
        self.sensor_type, = struct.unpack(self.endian + 'B', dis.read(1))
        self.option_flag, = struct.unpack(self.endian + 'B', dis.read(1))
        self.site_name = dis.read(5).decode("utf-8").strip()
        # self.site_name = self.site_name.encode('ascii',errors='ignore').decode()
        self.site_name = ''.join([s for s in self.site_name if ord(s) < 127 and ord(s)>0])
        self.channel_name = dis.read(3).decode("utf-8").strip()
        self.location_name = dis.read(2).decode("utf-8").strip()
        self.data_format = dis.read(2).decode("utf-8").strip()
        self.calib_factor, = struct.unpack(self.endian + 'f', dis.read(4))
        self.calib_period, = struct.unpack(self.endian + 'f', dis.read(4))
        self.time_stamp = dis.read(20).decode("utf-8").strip()

        self.subframe_time_length, = struct.unpack(self.endian + 'I', dis.read(4))
        self.npts, = struct.unpack(self.endian + 'I', dis.read(4))
        self.channel_status_size, = struct.unpack(self.endian + 'I', dis.read(4))
        self.channel_status_data = dis.read(self.channel_status_size)
        self.data_size, = struct.unpack(self.endian + 'I', dis.read(4))
        pad_length = (4 - self.data_size % 4) % 4
        fmt = "{0}B".format(self.endian)
        self.channel_data = [struct.unpack(fmt, dis.read(1))[0] for p in range(self.data_size)]
        # self.channel_data = dis.read(self.data_size)
        dis.read(pad_length)
        self.subframe_count, = struct.unpack(self.endian + 'I', dis.read(4))
        self.auth_key_ident, = struct.unpack(self.endian + 'I', dis.read(4))
        self.auth_size, = struct.unpack(self.endian + 'I', dis.read(4))
        self.auth_value = dis.read(self.auth_size)

        self.next_value = np.zeros(1)
        self.data = uncompress(self.channel_data, self.npts, self.next_value, self.endian)
        self.next_sample = self.next_value[0]
        trase_header = Stats()
        trase_header.station = self.site_name
        trase_header.channel = self.channel_name
        trase_header.network = self.location_name
        trase_header.starttime = UTCDateTime.strptime(self.time_stamp, "%Y%j %H:%M:%S.%f")
        trase_header.sampling_rate = self.npts / (self.subframe_time_length/1000.0)
        trase_header.npts = self.npts
        self.trace = Trace(data=self.data, header=trase_header)
        log.debug(self.trace)
        # self.trase.plot()


class CD11Frame(object):
    endian = ">"
    frameType = 5
    trailer_offset = 0
    srcName = ""
    dstName = ""
    sequenceNumber = 0
    series = 0

    def __init__(self, order=">"):
        self.endian = order

    def read(self, dis):
        self.frameType, = struct.unpack(self.endian + 'I', dis.read(4))
        self.trailer_offset, = struct.unpack(self.endian + 'I', dis.read(4))
        self.srcName = dis.read(8).decode("utf-8").strip()
        self.dstName = dis.read(8).decode("utf-8").strip()
        self.sequenceNumber, = struct.unpack(self.endian + 'Q', dis.read(8))
        self.series, = struct.unpack(self.endian + 'I', dis.read(4))


class CD11(object):
    endian = ">"
    numChans = 0
    frameDurnMillis = 0
    padding = None
    channelString = None
    frameTimeStr = None
    nominalTimeMillis = None
    subFrames = []
    tx_state = None
    check_order = True
    buffer = None

    def __init__(self, order=">", check_order=True):
        self.endian = order
        self.check_order = check_order
        self.header = CD11Frame(self.endian)
        self.trailer = FrameTrailer(self.endian)

    def check_endian(self, buff):
        # dis = open(filename, "rb")
        # b = dis.read(4)
        # dis.close()
        type = struct.unpack('<I', buff)[0]
        if type == 5:
            log.debug("Type: BIG_ENDIAN")
            self.endian = "<"
        else:
            type = struct.unpack('>I', buff)[0]
            if type == 5:
                log.debug("Type: LITTLE_ENDIAN")
                self.endian = ">"
        return self.endian

    def read(self, filename):
        stream = Stream()
        dis = open(filename, "rb")
        buffer = dis.read()
        dis.close()
        if self.check_order:
            self.check_endian(buffer[:4])
        f = io.BytesIO(buffer)
        self.header.read(f)
        rez = self.read_frames(f, stream)
        self.trailer.read(f)
        f.close()
        stream.merge()
        return stream

    def read_frames(self, dis, stream):
        self.numChans, = struct.unpack(self.endian + 'I', dis.read(4))
        self.frameDurnMillis, = struct.unpack(self.endian + 'I', dis.read(4))
        self.frameTimeStr = dis.read(20).decode("utf-8").strip()
        self.nominalTimeMillis = UTCDateTime.strptime(self.frameTimeStr, "%Y%j %H:%M:%S.%f")
        ch_str_len, = struct.unpack(self.endian + 'I', dis.read(4))
        try:
            self.channelString = dis.read(ch_str_len).decode("utf-8").strip()
        except:
            pass
        pad_length = (4 - ch_str_len % 4) % 4
        fmt = "{0}B".format(self.endian)
        self.padding = [struct.unpack(fmt, dis.read(1))[0] for p in range(pad_length)]
        self.subFrames = []
        try:
            while True:
                chn_s_frame = ChannelSubFrame(order=self.endian)
                chn_s_frame.read(dis)
                self.subFrames.append(chn_s_frame)
                stream.append(chn_s_frame.trace)
        except:
            return False
        return True
