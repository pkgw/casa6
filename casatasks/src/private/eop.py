import math
import re
import shutil

from casatasks import casalog
from casatools import ctsys
from casatools import table
from casatools import msmetadata

import numpy as np
import erfa

c = 299792458.0
arcsec = 0.00000484813681109535993589914098765

def ut1_utc(eops, jd1, jd2):
    mjd1, mjd2 = jd1 - 2400000.5, jd2
    i = np.searchsorted(eops['MJD'], mjd1 + mjd2, side='right')
    i1 = np.clip(i, 1, len(eops) - 1)
    i0 = i1 - 1
    mjd_0, mjd_1 = eops['MJD'][i0], eops['MJD'][i1]
    val_0, val_1 = eops['UT1_UTC'][i0], eops['UT1_UTC'][i1]
    d_val = val_1 - val_0
    d_val -= round(d_val)
    val = val_0 + (mjd1 - mjd_0 + mjd2) / (mjd_1 - mjd_0) * d_val
    return val

def pm_xy(eops, jd1, jd2):
    mjd1, mjd2 = jd1 - 2400000.5, jd2
    i = np.searchsorted(eops['MJD'], mjd1 + mjd2, side='right')
    i1 = np.clip(i, 1, len(eops) - 1)
    i0 = i1 - 1
    mjd_0, mjd_1 = eops['MJD'][i0], eops['MJD'][i1]
    val_0, val_1 = eops['PM_x'][i0], eops['PM_x'][i1]
    xp = val_0 + (mjd1 - mjd_0 + mjd2) / (mjd_1 - mjd_0) * (val_1 - val_0)
    val_0, val_1 = eops['PM_y'][i0], eops['PM_y'][i1]
    yp = val_0 + (mjd1 - mjd_0 + mjd2) / (mjd_1 - mjd_0) * (val_1 - val_0)
    return xp * arcsec, yp * arcsec

def utctt(jd1, jd2):
    tai = erfa.utctai(jd1, jd2)
    return erfa.taitt(tai[0], tai[1])

def calc_delay(eops, pos, p, mjd):
    xp, yp = pm_xy(eops, 2400000.5, mjd)
    dut1 = ut1_utc(eops, 2400000.5, mjd) / 86400
    tt = utctt(2400000.5, mjd)
    r = erfa.c2t00a(tt[0], tt[1], 2400000.5 + mjd, dut1, xp, yp)
    p = erfa.rxp(r, p)
    return erfa.pdp(pos, p) / c

def add_row(tb, obsid, field, antenna, spws, reffreq, t, pos, pc,
            old_eops, new_eops, delta):
    dtm = delta / 86400
    tm = t / 86400
    old_delay = [0, 0, 0]
    old_delay[0] = calc_delay(old_eops, pos, pc, tm - dtm)
    old_delay[1] = calc_delay(old_eops, pos, pc, tm)
    old_delay[2] = calc_delay(old_eops, pos, pc, tm + dtm)
    new_delay = [0, 0, 0]
    new_delay[0] = calc_delay(new_eops, pos, pc, tm - dtm)
    new_delay[1] = calc_delay(new_eops, pos, pc, tm)
    new_delay[2] = calc_delay(new_eops, pos, pc, tm + dtm)
    delay = new_delay[1] - old_delay[1]
    new_rate = (new_delay[2] - new_delay[0]) / (2 * delta)
    old_rate = (old_delay[2] - old_delay[0]) / (2 * delta)
    rate = new_rate - old_rate

    param = np.zeros(shape=(8,1), dtype='float32')
    paramerr = -np.ones(shape=(8,1), dtype='float32')
    flag = np.zeros(shape=(8,1), dtype='float32')
    snr = np.ones(shape=(8,1), dtype='float32')
    weight = np.ones(shape=(8,1), dtype='float32')

    param[1] = delay * 1e9
    param[2] = rate
    param[5] = delay * 1e9
    param[6] = rate

    for spw in spws:
        param[0] = 2 * math.pi * delay * reffreq[spw]
        param[4] = 2 * math.pi * delay * reffreq[spw]

        row = tb.nrows()
        tb.addrows(1)
        tb.putcell('TIME', row, t)
        tb.putcell('INTERVAL', row, 0)
        tb.putcell('ANTENNA1', row, antenna)
        tb.putcell('ANTENNA2', row, -1)
        tb.putcell('FIELD_ID', row, field)
        tb.putcell('SCAN_NUMBER', row, -1)
        tb.putcell('OBSERVATION_ID', row, obsid)
        tb.putcell('SPECTRAL_WINDOW_ID', row, spw)
        tb.putcell('FPARAM', row, param)
        tb.putcell('PARAMERR', row, paramerr)
        tb.putcell('FLAG', row, flag)
        tb.putcell('SNR', row, snr)
        tb.putcell('WEIGHT', row, weight)
        continue
    pass

def get_eops_from_ms(vis, obsid):
    tb = table()

    eops = {}
    eops['MJD'] = []
    eops['UT1_UTC'] = []
    eops['PM_x'] = []
    eops['PM_y'] = []

    try:
        tb.open(vis + '/EARTH_ORIENTATION')
        res = tb.query('OBSERVATION_ID == %d' % obsid)
        for row in range(res.nrows()):
            pm = res.getcell('PM', row)
            eops['MJD'].append(res.getcell('TIME', row) / 86400)
            eops['PM_x'].append(pm[0] / arcsec)
            eops['PM_y'].append(pm[1] / arcsec)
            eops['UT1_UTC'].append(res.getcell('UT1_UTC', row))
            continue
        tb.close()
    except:
        return None

    if len(eops['MJD']) == 0:
        return None
    return eops

def get_eops_from_casadata(mjd_min, mjd_max):
    tb = table()

    eops = {}
    eops['MJD'] = []
    eops['UT1_UTC'] = []
    eops['PM_x'] = []
    eops['PM_y'] = []

    try:
        path = ctsys.resolve('geodetic/IERSeop2000')
        tb.open(path)
        res = tb.query('MJD > %f && MJD < %f' % (mjd_min - 0.5, mjd_max + 0.5))
        for row in range(res.nrows()):
            eops['MJD'].append(res.getcell('MJD', row))
            eops['UT1_UTC'].append(res.getcell('dUT1', row))
            eops['PM_x'].append(res.getcell('x', row))
            eops['PM_y'].append(res.getcell('y', row))
            continue
        tb.close()
    except:
        return None

    if len(eops['MJD']) == 0:
        return None
    return eops

# Parse UNSO finals version 2.0 format
def parse_usno_2_0(fp, mjd_min, mjd_max):
    eops = {}
    eops['MJD'] = []
    eops['UT1_UTC'] = []
    eops['PM_x'] = []
    eops['PM_y'] = []

    for line in fp:
        if re.match(r'#', line):
            continue
        jd_tai = float(line[0:9])
        mjd = jd_tai - 2400000.5
        if mjd < mjd_min:
            continue
        if mjd > mjd_max:
            continue
        pm_x = float(line[10:17]) / 10
        pm_y = float(line[18:25]) / 10
        ut1_tai = float(line[26:35]) / 1e6
        jd_utc = erfa.taiutc(jd_tai, 0)
        utc_tai = (jd_utc[1] + (jd_utc[0] - jd_tai)) * 86400
        ut1_utc =  ut1_tai - utc_tai
        mjd = jd_utc[1] + (jd_utc[0] - 2400000.5)
        eops['MJD'].append(mjd)
        eops['UT1_UTC'].append(ut1_utc)
        eops['PM_x'].append(pm_x)
        eops['PM_y'].append(pm_y)
        continue

    if len(eops['MJD']) == 0:
        return None
    return eops

# Parse UNSO finals version 2.1 format
def parse_usno_2_1(fp, mjd_min, mjd_max):
    eops = {}
    eops['MJD'] = []
    eops['UT1_UTC'] = []
    eops['PM_x'] = []
    eops['PM_y'] = []

    for line in fp:
        if re.match(r'#', line):
            continue
        jd_tai = float(line[0:9])
        mjd = jd_tai - 2400000.5
        if mjd < mjd_min:
            continue
        if mjd > mjd_max:
            continue
        pm_x = float(line[10:18]) / 10
        pm_y = float(line[19:27]) / 10
        ut1_tai = float(line[28:39]) / 1e6
        jd_utc = erfa.taiutc(jd_tai, 0)
        utc_tai = (jd_utc[1] + (jd_utc[0] - jd_tai)) * 86400
        ut1_utc =  ut1_tai - utc_tai
        mjd = jd_utc[1] + (jd_utc[0] - 2400000.5)
        eops['MJD'].append(mjd)
        eops['UT1_UTC'].append(ut1_utc)
        eops['PM_x'].append(pm_x)
        eops['PM_y'].append(pm_y)
        continue

    if len(eops['MJD']) == 0:
        return None
    return eops

# Parse old IERS format
def parse_eopc04(fp, mjd_min, mjd_max):
    eops = {}
    eops['MJD'] = []
    eops['UT1_UTC'] = []
    eops['PM_x'] = []
    eops['PM_y'] = []

    for line in fp:
        if re.match(r'\s', line) or re.match(r'#', line):
            continue
        mjd = int(line[12:19])
        if mjd < mjd_min:
            continue
        if mjd > mjd_max:
            continue
        pm_x = float(line[20:30])
        pm_y = float(line[31:41])
        ut1_utc = float(line[42:53])
        eops['MJD'].append(mjd)
        eops['UT1_UTC'].append(ut1_utc)
        eops['PM_x'].append(pm_x)
        eops['PM_y'].append(pm_y)
        continue

    if len(eops['MJD']) == 0:
        return None
    return eops

def get_eops_from_file(infile, mjd_min, mjd_max):
    fp = open(infile)
    if not fp:
        msg = 'Cannot open ' + infile
        casalog.post(msg, 'SEVERE')
        return None
    try:
        line = fp.readline()
        if line.startswith("EOP-MOD Ver 2.0"):
            return parse_usno_2_0(fp, mjd_min, mjd_max)
        elif line.startswith("EOP-MOD Ver 2.1"):
            return parse_usno_2_1(fp, mjd_min, mjd_max)
        else:
            return parse_eopc04(fp, mjd_min, mjd_max)
    except:
        return None

def do_generate_eop(vis, caltable, infile):
    msmd = msmetadata()
    tb = table()

    # Use a 1 minute calibration interval.  This should be sufficient
    # as we'll use both the delay and delay rate in the interpolation.
    delta = 60

    tb.open(caltable + '/SPECTRAL_WINDOW')
    reffreq = tb.getcol('CHAN_FREQ')[0]
    tb.close()

    msmd.open(vis)
    tb.open(caltable, nomodify=False)

    # Prefetch antenna positions.
    positions = {}
    for antenna in msmd.antennaids():
        pos = msmd.antennaposition(antenna)
        pos = erfa.s2p(pos['m0']['value'], pos['m1']['value'], pos['m2']['value'])
        positions[antenna] = pos
        continue

    # Prefetch phase centers.
    phasecenters = {}
    for field in msmd.fieldsforname():
        pc = msmd.phasecenter(field)
        pc = erfa.s2c(pc['m0']['value'], pc['m1']['value'])
        phasecenters[field] = pc
        continue

    # The original EOPs may include differen values for the same day
    # for different observations.  Therefore we iterate over
    # observations and create calibration table entries that include
    # the OBSERVATION_ID.
    for obsid in range(msmd.nobservations()):
        # Get the original EOPs for this obervation from the MS.
        old_eops = get_eops_from_ms(vis, obsid)
        if not old_eops:
            msg = 'No EOPS found in ' + vis
            casalog.post(msg, 'SEVERE')
            continue

        # Get the updated EOPs for the timerange covered by the
        # original EOPs.
        mjd_min = np.min(old_eops['MJD'])
        mjd_max = np.max(old_eops['MJD'])
        if infile and infile != 'None':
            new_eops = get_eops_from_file(infile, mjd_min, mjd_max)
        else:
            new_eops = get_eops_from_casadata(mjd_min, mjd_max)
        if not new_eops:
            msg = 'No updated EOPS found for MJD ' + str(mjd_min) + '-' \
                + str(mjd_max) + ' in ' + infile
            casalog.post(msg, 'SEVERE')
            continue

        # Iterate by scan over all the antennas and fields in this
        # observation.
        for scan in msmd.scannumbers(obsid=obsid):
            for antenna in msmd.antennasforscan(scan=scan, obsid=obsid):
                for field in msmd.fieldsforscan(scan=scan, obsid=obsid):
                    pos = positions[antenna]
                    pc = phasecenters[field]

                    spws = msmd.spwsforscan(scan=scan, obsid=obsid)
                    times = msmd.timesforscan(scan=scan, obsid=obsid)
                    t0 = 0
                    for t in times:
                        if t - t0 > delta:
                            t0 = t
                            add_row(tb, obsid, field, antenna, spws, reffreq,
                                    t, pos, pc, old_eops, new_eops, delta)
                            pass
                        continue
                    if t != t0:
                        add_row(tb, obsid, field, antenna, spws, reffreq,
                                t, pos, pc, old_eops, new_eops, delta)
                        pass

                    continue
                continue
            continue
        continue

    nrows = tb.nrows()
    tb.close()

    return nrows

def generate_eop(vis, caltable, infile):
    nrows = 0
    try:
        nrows = do_generate_eop(vis, caltable, infile)
    finally:
        if nrows == 0:
            shutil.rmtree(caltable)
