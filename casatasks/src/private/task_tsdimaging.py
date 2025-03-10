# sd task for imaging

import collections
import contextlib
import os
import re
import shutil
import time

import numpy

from functools import partial

from casatasks import casalog
from casatools import image, imager
from casatools import ms as mstool
from casatools import quanta

from . import mslisthelper, sdbeamutil, sdutil
# (1) Import the python application layer
from .imagerhelpers.imager_base import PySynthesisImager
from .imagerhelpers.input_parameters import ImagerParameters

image_suffix = '.image'
residual_suffix = '.residual'
weight_suffix = '.weight'
associate_suffixes = ['.psf', '.sumwt', weight_suffix, residual_suffix]


@contextlib.contextmanager
def open_ia(imagename):
    ia = image()
    ia.open(imagename)
    try:
        yield ia
    finally:
        ia.close()


@contextlib.contextmanager
def open_ms(vis):
    ms = mstool()
    ms.open(vis)
    try:
        yield ms
    finally:
        ms.close()


class SelectionHandler(object):
    def __init__(self, sel):
        self.sel = sel
        if isinstance(self.sel, str):
            self.selector = self._select0
        elif len(self.sel) == 1:
            self.selector = self._select1
        else:
            self.selector = self._select2

    def __call__(self, i):
        return self.selector(i)

    def _select0(self, i):
        return self.sel

    def _select1(self, i):
        return self.sel[0]

    def _select2(self, i):
        return self.sel[i]


class OldImagerBasedTools(object):
    def __init__(self):
        self.imager = imager()

    @contextlib.contextmanager
    def open_old_imager(self, vis):
        try:
            self.imager.open(vis)
            yield self.imager
        finally:
            self.imager.close()

    @contextlib.contextmanager
    def open_and_select_old_imager(self, vislist, field, spw, antenna, scan, intent, timerange):
        if isinstance(vislist, str):
            with self.open_old_imager(vislist) as im:
                im.selectvis(field=field,
                             spw=spw,
                             nchan=-1,
                             start=0,
                             step=1,
                             baseline=antenna,
                             scan=scan,
                             intent=intent,
                             time=timerange)
                yield im
        else:
            fieldsel = SelectionHandler(field)
            spwsel = SelectionHandler(spw)
            antennasel = SelectionHandler(antenna)
            scansel = SelectionHandler(scan)
            intentsel = SelectionHandler(intent)
            timerangesel = SelectionHandler(timerange)
            try:
                for i in range(len(vislist)):
                    vis = vislist[i]
                    _field = fieldsel(i)
                    _spw = spwsel(i)
                    _antenna = antennasel(i)
                    _scan = scansel(i)
                    _intent = intentsel(i)
                    _timerangesel = timerangesel(i)
                    if len(_antenna) == 0:
                        _baseline = _antenna
                    elif len(_antenna) < 4 or _antenna[:-3] != '&&&':
                        _baseline = _antenna + '&&&'
                    else:
                        _baseline = _antenna
                    self.imager.selectvis(
                        vis, field=_field, spw=_spw, nchan=-1, start=0, step=1,
                        baseline=_baseline, scan=_scan, intent=_intent, time=_timerangesel)
                yield self.imager
            finally:
                self.imager.close()

    def test(self, vis):
        with self.open_old_imager(vis):
            casalog.post('test')
            raise RuntimeError('ERROR!')

    def get_pointing_sampling_params(self, vis, field, spw, baseline, scan, intent, timerange,
                                     outref, movingsource, pointingcolumntouse, antenna_name):
        with self.open_old_imager(vis) as im:
            im.selectvis(field=field,
                         spw=spw,
                         nchan=-1,
                         start=0,
                         step=1,
                         baseline=baseline,
                         scan=scan,
                         intent=intent,
                         time=timerange)
            sampling_params = im.pointingsampling(pattern='raster',
                                                  ref=outref,
                                                  movingsource=movingsource,
                                                  pointingcolumntouse=pointingcolumntouse,
                                                  antenna='{0}&&&'.format(antenna_name))
        return sampling_params

    def get_map_extent(self, vislist, field, spw, antenna, scan, intent, timerange,
                       ref, movingsource, pointingcolumntouse):

        with self.open_and_select_old_imager(vislist=vislist, field=field,
                                             spw=spw, antenna=antenna, scan=scan,
                                             intent=intent, timerange=timerange) as im:
            map_param = im.mapextent(ref=ref, movingsource=movingsource,
                                     pointingcolumntouse=pointingcolumntouse)
        return map_param


def check_conformance(mslist, check_result):
    """Check conformance of input MS list

    Check conformance of input MS list, particularlly existence of
    WEIGHT_SPECTRUM column.

    Args:
        mslist (list): list of names for input MS
        check_result (dict): result of conformance check.
            see mslisthelper.check_mslist for detail about
            the structure of check_result.

    Returns:
        dict: Per-column set of names for MS that needs to be
              edited to resolve the conformance. Top level dict
              has two keys, "remove" and "add", which indicate
              the operation to be applied to the columns.
    """
    process_dict = {
        'remove': collections.defaultdict(set),
        'add': collections.defaultdict(set)
    }
    for name, summary in check_result.items():
        if 'Main' in summary:
            missingcol_list = [summary['Main']['missingcol_{}'.format(x)] for x in ['a', 'b']]
            ms_a = mslist[0]
            ms_b = name

            # "remove" operation:
            #  - MS list is opposite order to missingcol_list
            #  - if WEIGHT_SPECTRUM column is missing in one MS,
            #    the column should be removed from another MS.
            for c, t in zip(missingcol_list, [ms_b, ms_a]):
                if 'WEIGHT_SPECTRUM' in c:
                    process_dict['remove']['WEIGHT_SPECTRUM'].add(t)
            # "add" operation:
            #  - MS list is same order as missingcol_list
            #  - if CORRECTED_DATA column is missing in one MS,
            #    the column should be added to that MS.
            for c, t in zip(missingcol_list, [ms_a, ms_b]):
                if 'CORRECTED_DATA' in c:
                    process_dict['add']['CORRECTED_DATA'].add(t)
    return process_dict


def report_conformance(mslist, column_name, process_set):
    """Report conformance of input MS

    Report conformance of input MS, particularlly on the existence
    of the column given by column_name.

    Args:
        mslist (list): list of names for input MS
        column_name (str): name of the column
        process_set (set): set of names of MS that need to be edited
    """
    if len(process_set) > 0:
        casalog.post('Detected non-conformance of {} column in input list of MSes.'.format(column_name), priority='WARN')
        cols = ['exists?', 'MS name']
        header = ' '.join(cols)
        casalog.post('', priority='WARN')
        casalog.post('Summary of existence of {}:'.format(column_name), priority='WARN')
        casalog.post(header, priority='WARN')
        casalog.post('-' * len(header), priority='WARN')
        for name in mslist:
            basename = os.path.basename(name.rstrip('/'))
            exists = 'YES' if name in process_set else 'NO'
            row = '{:^7s} {:<s}'.format(exists, basename)
            casalog.post(row, priority='WARN')


def fix_conformance(process_dict):
    """Resolve non-conformance by removing WEIGHT_SPECTRUM

    Two non-conformances are fixed, WEIGHT_SPECTRUM and
    CORRECTED_DATA. Remove WEIGHT_SPECTRUM column from, or
    add CORRECTED_DATA to the MS provided by process_dict.
    Backup is created with the name:

      <original_name>.sdimaging.backup-<timestamp>

    Args:
        process_dict (dict): per-operation ("remove" and "add")
                             key-value pair of column name and
                             list of names for MS to be edited

    Returns:
        dict: mapping of original MS name and the name of backup
    """
    backup_list = {}
    process_list = set()
    for v in process_dict.values():
        for w in v.values():
            process_list = process_list.union(w)
    for name in process_list:
        basename = os.path.basename(name.rstrip('/'))
        timestamp = time.strftime('%Y%m%dT%H%M%S', time.gmtime())
        backup_name = basename + '.sdimaging.backup-{}'.format(timestamp)
        with sdutil.table_manager(name) as tb:
            tb.copy(backup_name, deep=True, returnobject=True).close()
        backup_list[name] = backup_name
        casalog.post('Copy of "{}" has been saved to "{}"'.format(name, backup_name), priority='WARN')

    for colname, msnames in process_dict['remove'].items():
        for name in msnames:
            casalog.post('{} will be removed from "{}"'.format(colname, name), priority='WARN')
            with sdutil.table_manager(name, nomodify=False) as tb:
                if colname in tb.colnames():
                    tb.removecols(colname)

    for colname, msnames in process_dict['add'].items():
        for name in msnames:
            casalog.post('{} will be added to "{}"'.format(colname, name), priority='WARN')
            with sdutil.calibrater_manager(name, addmodel=False, addcorr=True):
                pass
    return backup_list


def conform_mslist(mslist, ignore_columns=['CORRECTED_DATA']):
    """Make given set of MS data conform

    Here, only conformance on the existence of WEIGHT_SPECTRUM
    is checked and resolved because non-conformance of WEIGHT_SPECTRUM,
    i.e. some MS have the column while others don't, could cause
    the task to crash. If non-conformance is detected, all existing
    WEIGHT_SPECTRUM columns are removed. This opration modifies input
    MS so data will be backed up with the name:

      <original_name>.sdimaging.backup-<timestamp>

    Args:
        mslist (list): list of names for input MS
    """
    check_result = mslisthelper.check_mslist(mslist, testcontent=False)
    process_dict = check_conformance(mslist, check_result)
    fix_dict = dict()
    for op, op_dict in process_dict.items():
        fix_dict[op] = dict()
        for col, process_set in op_dict.items():
            if col not in ignore_columns:
                report_conformance(mslist, col, process_set)
                fix_dict[op][col] = process_set
    fix_conformance(fix_dict)


def sort_vis(vislist, spw, mode, width, field, antenna, scan, intent, timerange):
    """Sort the given MeasurementSet path(s) by their earliest data-taking time, in increasing order.

    Return a 7-tuple where
        * the first entry holds the re-ordered paths
        * the remaining entries hold their corresponding data selection parameters
        * FIXME: input parameters mode and width are not used
    """
    if isinstance(vislist, str) or len(vislist) == 1:
        return vislist, field, spw, antenna, scan, intent, timerange
    # chronological sort
    sorted_vislist, sorted_timelist = mslisthelper.sort_mslist(vislist)
    _vislist = list(vislist)
    sorted_idx = [_vislist.index(vis) for vis in sorted_vislist]
    mslisthelper.report_sort_result(sorted_vislist, sorted_timelist, sorted_idx, mycasalog=casalog)
    # conform MS
    conform_mslist(sorted_vislist)
    fieldsel = SelectionHandler(field)
    sorted_field = [fieldsel(i) for i in sorted_idx]
    spwsel = SelectionHandler(spw)
    sorted_spw = [spwsel(i) for i in sorted_idx]
    antennasel = SelectionHandler(antenna)
    sorted_antenna = [antennasel(i) for i in sorted_idx]
    scansel = SelectionHandler(scan)
    sorted_scan = [scansel(i) for i in sorted_idx]
    intentsel = SelectionHandler(intent)
    sorted_intent = [intentsel(i) for i in sorted_idx]
    timerangesel = SelectionHandler(timerange)
    sorted_timerange = [timerangesel(i) for i in sorted_idx]
    return sorted_vislist, sorted_field, sorted_spw, sorted_antenna, sorted_scan, sorted_intent, sorted_timerange


def _configure_spectral_axis(mode, nchan, start, width, restfreq):
    # fix default
    if mode == 'channel':
        if start == '':
            start = 0
        if width == '':
            width = 1
    else:
        if start == 0:
            start = ''
        if width == 1:
            width = ''
    # fix unit
    if mode == 'frequency':
        myunit = 'Hz'
    elif mode == 'velocity':
        myunit = 'km/s'
    else:  # channel
        myunit = ''

    tmp_start = _format_quantum_unit(start, myunit)
    if tmp_start is None:
        raise ValueError("Invalid unit for %s in mode %s: %s" % ('start', mode, start))
    start = tmp_start
    if mode == 'channel':
        start = int(start)
    tmp_width = _format_quantum_unit(width, myunit)
    if tmp_width is None:
        raise ValueError("Invalid unit for %s in mode %s: %s" % ('width', mode, width))
    width = tmp_width
    if mode == 'channel':
        width = int(width)

    # TODO: work for nchan
    imnchan = nchan
    imstart = start
    imwidth = width
    return imnchan, imstart, imwidth


def _format_quantum_unit(data, unit):
    """Format quantity data.

    Returns False if data has an unit which in not a variation of
    input unit.
    Otherwise, returns input data as a quantum string. The input
    unit is added to the return value if no unit is in data.
    """
    my_qa = quanta()
    if data == '' or my_qa.compare(data, unit):
        return data
    if my_qa.getunit(data) == '':
        casalog.post("No unit specified. Using '%s'" % unit)
        return '%f%s' % (data, unit)
    return None


def _handle_grid_defaults(value):
    ret = ''
    if isinstance(value, int) or isinstance(value, float):
        ret = str(value)
    elif isinstance(value, str):
        ret = value
    return ret


def _calc_PB(vis, antenna_id, restfreq):
    """Calculate the primary beam size of antenna.

    Calculate the primary beam size of antenna, using dish diamenter
    and rest frequency
    Average antenna diamter and reference frequency are adopted for
    calculation.
    The input argument should be a list of antenna IDs.
    """
    logger = sdutil.Casalog(origin="_calc_PB")
    logger.post("Calculating Primary beam size:")
    # CAS-5410 Use private tools inside task scripts
    my_qa = quanta()

    pb_factor = 1.175
    # Reference frequency
    ref_freq = restfreq
    if type(ref_freq) in [float, numpy.float64]:
        ref_freq = my_qa.tos(my_qa.quantity(ref_freq, 'Hz'))
    if not my_qa.compare(ref_freq, 'Hz'):
        msg = "Could not get the reference frequency. " + \
              "Your data does not seem to have valid one in selected field.\n" + \
              "PB is not calculated.\n" + \
              "Please set restreq or cell manually to generate an image."
        raise RuntimeError(msg)
    # Antenna diameter
    with sdutil.table_manager(os.path.join(vis, 'ANTENNA')) as tb:
        antdiam_ave = tb.getcell('DISH_DIAMETER', antenna_id)
    # antdiam_ave = self._get_average_antenna_diameter(antenna)
    # Calculate PB
    wave_length = 0.2997924 / my_qa.convert(my_qa.quantity(ref_freq), 'GHz')['value']
    D_m = my_qa.convert(antdiam_ave, 'm')['value']
    lambda_D = wave_length / D_m * 3600. * 180 / numpy.pi
    PB = my_qa.quantity(pb_factor * lambda_D, 'arcsec')
    # Summary
    logger.post(f"- Antenna diameter: {D_m} m")
    logger.post(f"- Reference Frequency: {ref_freq}")
    logger.post(f"PB size = {pb_factor:5.3f} * lambda/D = {my_qa.tos(PB)}")
    return PB


def _get_imsize(width, height, dx, dy):
    casalog.post("Calculating pixel size.")
    # CAS-5410 Use private tools inside task scripts
    my_qa = quanta()
    ny = numpy.ceil((my_qa.convert(height, my_qa.getunit(dy))['value'] /
                     my_qa.getvalue(dy)))
    nx = numpy.ceil((my_qa.convert(width, my_qa.getunit(dx))['value'] /
                     my_qa.getvalue(dx)))
    casalog.post("- Map extent: [%s, %s]" % (my_qa.tos(width), my_qa.tos(height)))
    casalog.post("- Cell size: [%s, %s]" % (my_qa.tos(dx), my_qa.tos(dy)))
    casalog.post("Image pixel numbers to cover the extent: [%d, %d] (projected)" %
                 (nx + 1, ny + 1))
    return [int(nx + 1), int(ny + 1)]


def _get_pointing_extent(phasecenter, vislist, field, spw, antenna, scan, intent, timerange,
                         pointingcolumntouse, ephemsrcname):
    # MS selection is ignored. This is not quite right.
    casalog.post("Calculating map extent from pointings.")
    # CAS-5410 Use private tools inside task scripts
    my_qa = quanta()
    ret_dict = {}

    if isinstance(vislist, str):
        vis = vislist
    else:
        vis = vislist[0]

    # colname = pointingcolumntouse.upper()

    if phasecenter == "":
        # defaut is J2000
        base_mref = 'J2000'
    elif isinstance(phasecenter, int) or phasecenter.isdigit():
        # may be field id
        with sdutil.table_manager(os.path.join(vis, 'FIELD')) as tb:
            base_mref = tb.getcolkeyword('PHASE_DIR', 'MEASINFO')['Ref']
    else:
        # may be phasecenter is explicitly specified
        # numeric value: 3.14, -.3e1, etc.
        numeric_pattern = r'[-+]?([0-9]+(.[0-9]*)?|\.[0-9]+)([eE]-?[0-9])?'
        # HMS string: 9:15:29, -9h15m29
        hms_pattern = r'[-+]?[0-9]+[:h][0-9]+[:m][0-9.]+s?'
        # DMS string: 9.15.29, -9d15m29s
        dms_pattern = r'[-+]?[0-9]+[.d][0-9]+[.m][0-9.]+s?'
        # composite pattern
        pattern = fr'^({numeric_pattern}|{hms_pattern}|{dms_pattern})$'
        items = phasecenter.split()
        base_mref = 'J2000'
        for i in items:
            s = i.strip()
            if re.match(pattern, s) is None:
                base_mref = s
                break

    t = OldImagerBasedTools()
    mapextent = t.get_map_extent(vislist, field, spw, antenna, scan, intent, timerange,
                                 ref=base_mref, movingsource=ephemsrcname,
                                 pointingcolumntouse=pointingcolumntouse)
    # mapextent = self.imager.mapextent(ref=base_mref, movingsource=ephemsrcname,
    #                                  pointingcolumntouse=colname)
    if mapextent['status']:
        qheight = my_qa.quantity(mapextent['extent'][1], 'rad')
        qwidth = my_qa.quantity(mapextent['extent'][0], 'rad')
        qcent0 = my_qa.quantity(mapextent['center'][0], 'rad')
        qcent1 = my_qa.quantity(mapextent['center'][1], 'rad')
        scenter = '%s %s %s' % (base_mref, my_qa.formxxx(qcent0, 'hms'),
                                my_qa.formxxx(qcent1, 'dms'))

        casalog.post("- Pointing center: %s" % scenter)
        casalog.post("- Pointing extent: [%s, %s] (projected)" % (my_qa.tos(qwidth),
                                                                  my_qa.tos(qheight)))
        ret_dict['center'] = scenter
        ret_dict['width'] = qwidth
        ret_dict['height'] = qheight
    else:
        casalog.post(
            'Failed to derive map extent from the MSs registered to the imager probably '
            'due to mising valid data.',
            priority='SEVERE')
        ret_dict['center'] = ''
        ret_dict['width'] = my_qa.quantity(0.0, 'rad')
        ret_dict['height'] = my_qa.quantity(0.0, 'rad')
    return ret_dict


def _handle_image_params(imsize, cell, phasecenter,
                         vislist, field, spw, antenna, scan, intent, timerange,
                         restfreq, pointingcolumntouse, ephemsrcname):
    logger = sdutil.Casalog(origin="_handle_image_params")
    # round-up imsize
    _imsize = sdutil.to_list(imsize, int) or sdutil.to_list(imsize, numpy.integer)
    if _imsize is None:
        _imsize = imsize if hasattr(imsize, '__iter__') else [imsize]
        _imsize = [int(numpy.ceil(v)) for v in _imsize]
        logger.post(
            "imsize is not integers. force converting to integer pixel numbers.",
            priority="WARN")
        logger.post("rounded-up imsize: %s --> %s" % (str(imsize), str(_imsize)))

    # calculate cell based on PB if it is not given
    _cell = cell
    if _cell == '' or _cell[0] == '':
        # calc PB
        if isinstance(vislist, str):
            vis = vislist
        else:
            vis = vislist[0]
        if isinstance(antenna, str):
            antsel = antenna
        else:
            antsel = antenna[0]
        if antsel == '':
            antenna_id = 0
        else:
            if len(antsel) > 3 and antsel[:-3] == '&&&':
                baseline = antsel
            else:
                baseline = antsel + '&&&'
            with open_ms(vis) as ms:
                ms.msselect({'baseline': baseline})
                ndx = ms.msselectedindices()
            antenna_id = ndx['antenna1'][0]
        grid_factor = 3.
        logger.post("The cell size will be calculated using PB size of antennas in the first MS")
        qpb = _calc_PB(vis, antenna_id, restfreq)
        _cell = '%f%s' % (qpb['value'] / grid_factor, qpb['unit'])
        logger.post("Using cell size = PB/%4.2F = %s" % (grid_factor, _cell))

    # Calculate Pointing center and extent (if necessary)
    _phasecenter = phasecenter
    if _phasecenter == '' or len(_imsize) == 0 or _imsize[0] < 1:
        # return a dictionary with keys 'center', 'width', 'height'
        map_param = _get_pointing_extent(_phasecenter, vislist, field, spw, antenna, scan, intent,
                                         timerange, pointingcolumntouse, ephemsrcname)
        # imsize
        (cellx, celly) = sdutil.get_cellx_celly(_cell, unit='arcmin')
        if len(_imsize) == 0 or _imsize[0] < 1:
            _imsize = _get_imsize(map_param['width'], map_param['height'], cellx, celly)
            if _phasecenter != "":
                logger.post(
                    "You defined phasecenter but not imsize. "
                    "The image will cover as wide area as pointing in MS extends, "
                    "but be centered at phasecenter. "
                    "This could result in a strange image if your phasecenter is "
                    "apart from the center of pointings",
                    priority='WARN')
            if _imsize[0] > 1024 or _imsize[1] > 1024:
                logger.post(
                    "The calculated image pixel number is larger than 1024. "
                    "It could take time to generate the image depending on your computer resource. "
                    "Please wait...",
                    priority='WARN')

        # phasecenter
        # if empty, it should be determined here...
        if _phasecenter == "":
            _phasecenter = map_param['center']

    return _imsize, _cell, _phasecenter


def _get_param(ms_index, param):
    if isinstance(param, str):
        return param
    elif hasattr(param, '__iter__'):
        if len(param) == 1:
            return param[0]
        else:
            return param[ms_index]
    else:
        raise RuntimeError('Invalid parameter')


def _remove_image(imagename):
    if os.path.exists(imagename):
        if os.path.isdir(imagename):
            shutil.rmtree(imagename)
        elif os.path.isfile(imagename):
            os.remove(imagename)
        else:
            # could be a symlink
            os.remove(imagename)


def _get_restfreq_if_empty(vislist, spw, field, restfreq):
    qa = quanta()
    rf = None
    # if restfreq is nonzero float value, return it
    if isinstance(restfreq, float):
        if restfreq != 0.0:
            rf = restfreq
    # if restfreq is valid frequency string, return it
    # numeric string is interpreted as a value in the unit of Hz
    elif isinstance(restfreq, str):
        q = qa.convert(qa.quantity(restfreq), 'Hz')
        if q['unit'] == 'Hz' and q['value'] > 0.0:
            rf = restfreq
    # if restfreq is valid quantity, return it
    elif isinstance(restfreq, dict):
        q = qa.convert(restfreq, 'Hz')
        if q['unit'] == 'Hz' and q['value'] > 0.0:
            rf = restfreq

    if isinstance(vislist, str):
        vis = vislist
    elif hasattr(vislist, '__iter__'):
        vis = vislist[0]
    else:
        raise RuntimeError(
            'Internal Error: invalid vislist \'{0}\''.format(vislist))

    if isinstance(spw, str):
        spwsel = spw
    elif hasattr(spw, '__iter__'):
        spwsel = spw[0]
    else:
        raise RuntimeError(
            'Internal Error: invalid spw selection \'{0}\''.format(spw))

    if isinstance(field, str):
        fieldsel = field
    elif hasattr(field, '__iter__'):
        fieldsel = field[0]
    else:
        raise RuntimeError('Internal Error: invalid field selection \'{0}\''.format(field))

    with open_ms(vis) as ms:
        ms.msselect({'spw': spwsel, 'field': fieldsel})
        ndx = ms.msselectedindices()
        if len(ndx['spw']) > 0:
            spwid = ndx['spw'][0]
        else:
            spwid = None
        if len(ndx['field']) > 0:
            fieldid = ndx['field'][0]
        else:
            fieldid = None
    sourceid = None
    if fieldid is not None:
        with sdutil.table_manager(os.path.join(vis, 'FIELD')) as tb:
            sourceid = tb.getcell('SOURCE_ID', fieldid)
        if sourceid < 0:
            sourceid = None
    if rf is None:
        # if restfrequency is defined in SOURCE table, return it
        with sdutil.table_manager(os.path.join(vis, 'SOURCE')) as tb:
            if 'REST_FREQUENCY' in tb.colnames():
                tsel = None
                taql = ''
                if spwid is not None:
                    taql = 'SPECTRAL_WINDOW_ID == {0}'.format(spwid)
                if sourceid is not None:
                    delimiter = '&&' if len(taql) > 0 else ''
                    taql += '{0}SOURCE_ID == {1}'.format(delimiter, sourceid)
                if len(taql) > 0:
                    tsel = tb.query(taql)
                    t = tsel
                else:
                    t = tb
                try:
                    nrow = t.nrows()
                    if nrow > 0:
                        for irow in range(nrow):
                            if t.iscelldefined('REST_FREQUENCY', irow):
                                rfs = t.getcell('REST_FREQUENCY', irow)
                                if len(rfs) > 0:
                                    rf = rfs[0]
                                    break
                finally:
                    if tsel is not None:
                        tsel.close()

    if rf is None:
        if spwid is None:
            spwid = 0
        # otherwise, return mean frequency of given spectral window
        with sdutil.table_manager(os.path.join(vis, 'SPECTRAL_WINDOW')) as tb:
            cf = tb.getcell('CHAN_FREQ', spwid)
            rf = cf.mean()

    assert rf is not None

    return rf


def set_beam_size(vis, imagename,
                  field, spw, baseline, scan, intent, timerange,
                  ephemsrcname, pointingcolumntouse,
                  antenna_name, antenna_diameter,
                  restfreq,
                  gridfunction, convsupport, truncate, gwidth, jwidth):
    """Set estimated beam size to the image."""
    is_alma = antenna_name[0:2] in ['PM', 'DV', 'DA', 'CM']
    blockage = '0.75m' if is_alma else '0.0m'
    log_origin = 'set_beam_size'

    with open_ia(imagename) as ia:
        csys = ia.coordsys()
        outref = csys.referencecode('direction')[0]
        cell = list(csys.increment(type='direction', format='s')['string'])
        csys.done()

    old_tool = OldImagerBasedTools()
    sampling_params = old_tool.get_pointing_sampling_params(
        vis, field, spw, baseline,
        scan, intent, timerange,
        outref=outref,
        movingsource=ephemsrcname,
        pointingcolumntouse=pointingcolumntouse,
        antenna_name=antenna_name)
    qa = quanta()
    casalog.post(
        f'sampling_params={sampling_params}',
        origin=log_origin
    )
    xsampling, ysampling = qa.getvalue(qa.convert(sampling_params['sampling'],
                                       'arcsec'))
    angle = qa.getvalue(qa.convert(sampling_params['angle'], 'deg'))[0]

    casalog.post(
        f'Detected raster sampling = [{xsampling:f}, {ysampling:f}] arcsec',
        origin=log_origin
    )

    # handling of failed sampling detection
    valid_sampling = True
    # TODO: copy from sdimaging implementation
    sampling = [xsampling, ysampling]
    if abs(xsampling) < 2.2e-3 or not numpy.isfinite(xsampling):
        casalog.post(
            f"Invalid sampling={xsampling} arcsec. "
            f"Using the value of orthogonal direction={ysampling} arcsec",
            priority="WARN",
            origin=log_origin
        )
        sampling = [ysampling]
        angle = 0.0
        valid_sampling = False
    if abs(ysampling) < 1.0e-3 or not numpy.isfinite(ysampling):
        if valid_sampling:
            casalog.post(
                f"Invalid sampling={ysampling} arcsec. "
                f"Using the value of orthogonal direction={xsampling} arcsec",
                priority="WARN",
                origin=log_origin
            )
            sampling = [xsampling]
            angle = 0.0
            valid_sampling = True
    # reduce sampling and cell if it's possible
    if (len(sampling) > 1 and
            abs(sampling[0] - sampling[1]) <= 0.01 * abs(sampling[0])):
        sampling = [sampling[0]]
        angle = 0.0
        if cell[0] == cell[1]:
            cell = [cell[0]]
    if valid_sampling:
        # actual calculation of beam size
        bu = sdbeamutil.TheoreticalBeam()
        bu.set_antenna(antenna_diameter, blockage)
        bu.set_sampling(sampling, "%fdeg" % angle)
        bu.set_image_param(cell, restfreq, gridfunction,
                           convsupport, truncate, gwidth,
                           jwidth, is_alma)
        bu.summary()
        imbeam_dict = bu.get_beamsize_image()
        casalog.post(
            f"Setting image beam: "
            f"major={imbeam_dict['major']}, "
            f"minor={imbeam_dict['minor']}, "
            f"pa={imbeam_dict['pa']}",
            origin=log_origin
        )
        # set beam size to image
        with open_ia(imagename) as ia:
            ia.setrestoringbeam(**imbeam_dict)
    else:
        # BOTH sampling was invalid
        casalog.post(
            "Could not detect valid raster sampling. "
            "Exiting without setting beam size to image",
            priority='WARN',
            origin=log_origin
        )


def do_weight_mask(imagename, weightimage, minweight):
    # Mask image pixels whose weight are smaller than minweight.
    # Weight image should have 0 weight for pixels below < minweight
    logger = sdutil.Casalog(origin="do_weight_mask")
    logger.post(f"Start masking the map using minweight = {minweight:f}",
                priority="INFO")
    with open_ia(weightimage) as ia:
        try:
            stat = ia.statistics(mask="'" + weightimage + "' > 0.0",
                                 robust=True)
            valid_pixels = stat['npts']
        except RuntimeError as e:
            if 'No valid data found.' in str(e):
                valid_pixels = [0]
            else:
                raise e

    if len(valid_pixels) == 0 or valid_pixels[0] == 0:
        logger.post(
            "All pixels weight zero. "
            "This indicates no data in MS is in image area. "
            "Mask will not be set. Please check your image parameters.",
            priority="WARN")
        return
    median_weight = stat['median'][0]
    weight_threshold = median_weight * minweight
    logger.post(f"Median of weight in the map is {median_weight:f}",
                priority="INFO")
    logger.post("Pixels in map with weight <= median(weight)*minweight = "
                f"{weight_threshold:f} will be masked.",
                priority="INFO")
    # Leaving the original logic to calculate the number of masked pixels via
    # product of median of and min_weight (which i don't understand the logic)

    # Modify default mask
    with open_ia(imagename) as ia:
        ia.calcmask("'%s'>%f" % (weightimage, weight_threshold),
                    asdefault=True)

        ndim = len(ia.shape())
        _axes = numpy.arange(start=0 if ndim <= 2 else 2, stop=ndim)

        try:
            collapsed = ia.collapse('npts', axes=_axes)
            valid_pixels_after = collapsed.getchunk().sum()
            collapsed.close()
        except RuntimeError as e:
            if 'All selected pixels are masked' in str(e):
                valid_pixels_after = 0
            else:
                raise

    masked_fraction = 100. * (1. - valid_pixels_after / float(valid_pixels[0]))

    logger.post(f"This amounts to {masked_fraction:5.1f} % "
                "of the area with nonzero weight.",
                priority="INFO")
    logger.post(
        f"The weight image '{weightimage}' is returned by this task, "
        "if the user wishes to assess the results in detail.",
        priority="INFO")


def get_ms_column_unit(tb, colname):
    col_unit = ''
    if colname in tb.colnames():
        cdkw = tb.getcoldesc(colname)['keywords']
        for key in ['UNIT', 'QuantumUnits']:
            if key in cdkw:
                u = cdkw[key]
                if isinstance(u, str):
                    col_unit = u.strip()
                elif isinstance(u, (list, numpy.ndarray)) and len(u) > 0:
                    col_unit = u[0].strip()
            if col_unit:
                break
    return col_unit


def get_brightness_unit_from_ms(msname):
    image_unit = ''
    with sdutil.table_manager(msname) as tb:
        for column in ['CORRECTED_DATA', 'FLOAT_DATA', 'DATA']:
            image_unit = get_ms_column_unit(tb, column)
            if image_unit:
                break

    if image_unit.upper() == 'K':
        image_unit = 'K'
    else:
        image_unit = 'Jy/beam'

    return image_unit


@sdutil.sdtask_decorator
def tsdimaging(
        # Input data: list of MeasurementSets
        infiles,
        # Output data: CASA images: their path prefix, overwrite control
        outfile, overwrite,
        # Select data from input MeasurementSets, by
        field, spw, antenna, scan, intent, timerange,
        # Output images definition: frequency axis
        outframe,  # velocity frame
        mode, nchan, start, width, veltype,  # gridding type
        specmode,  # Doppler handling
        interpolation,  # interpolation mode
        # Output images definition: spatial axes
        pointingcolumn, convertfirst,
        projection,
        imsize, cell, phasecenter,
        # Output images definition: stokes axis
        stokes,
        # Gridder parameters
        gridfunction, convsupport, truncate, gwidth, jwidth,
        clipminmax,
        # Single-dish image: mask control
        minweight,
        # Single-dish image: metadata
        brightnessunit,
        # rest frequency to assign to image
        restfreq):

    origin = 'tsdimaging'
    imager = None

    try:  # Create the Single-Dish Image
        # Validate brightnessunit parameter CAS-11503
        image_unit = brightnessunit.lower().capitalize()
        if image_unit not in ['', 'K', 'Jy/beam']:
            raise ValueError(f"Invalid brightness unit: {brightnessunit}")

        # Handle outfile and overwrite parameters
        output_path_prefix = outfile.rstrip('/')
        singledish_image_path = output_path_prefix + image_suffix
        if os.path.exists(singledish_image_path):
            if not overwrite:
                raise RuntimeError(
                        f"Output file exists: '{singledish_image_path}'"
                      )
            else:
                # delete existing images
                if os.path.exists(singledish_image_path):
                    casalog.post(f"Removing '{singledish_image_path}'")
                    _remove_image(singledish_image_path)
                assert not os.path.exists(singledish_image_path)
                for _suffix in associate_suffixes:
                    path_to_remove = output_path_prefix + _suffix
                    if os.path.exists(path_to_remove):
                        casalog.post(f"Removing '{path_to_remove}'")
                        _remove_image(path_to_remove)
                    assert not os.path.exists(path_to_remove)

        # Tweak spw parameter into _spw
        if isinstance(spw, str):
            _spw = '*' + spw if spw.startswith(':') else spw
        else:
            _spw = ['*' + v if v.startswith(':') else v for v in spw]

        # Handle image spectral axis parameters
        imnchan, imstart, imwidth = _configure_spectral_axis(
            mode, nchan, start, width, restfreq
        )

        # Handle image restfreq parameter's default value
        _restfreq = _get_restfreq_if_empty(
            infiles, _spw, field, restfreq
        )

        # Handle gridder parameters
        # convert type of task's default values
        # to ones supported by the Synthesis Imager framework
        gtruncate = _handle_grid_defaults(truncate)
        ggwidth = _handle_grid_defaults(gwidth)
        gjwidth = _handle_grid_defaults(jwidth)

        # handle infiles parameter
        # sort input data to get results consistent with old sdimaging task
        _sorted = sort_vis(
            infiles, _spw, mode, imwidth, field,
            antenna, scan, intent, timerange
        )
        sorted_vis = _sorted[0]
        sorted_field = _sorted[1]
        sorted_spw = _sorted[2]
        sorted_antenna = _sorted[3]
        sorted_scan = _sorted[4]
        sorted_intent = _sorted[5]
        sorted_timerange = _sorted[6]

        def antenna_to_baseline(s):
            if len(s) == 0:
                return s
            elif len(s) > 3 and s.endswith('&&&'):
                return s
            else:
                return '{0}&&&'.format(s)

        if isinstance(sorted_antenna, str):
            sorted_baseline = antenna_to_baseline(sorted_antenna)
        else:
            sorted_baseline = [antenna_to_baseline(a) for a in sorted_antenna]

        # Handle image geometric parameters
        _ephemsrcname = ''
        ephem_sources = ['MERCURY', 'VENUS', 'MARS', 'JUPITER', 'SATURN',
                         'URANUS', 'NEPTUNE', 'PLUTO', 'SUN', 'MOON',
                         'TRACKFIELD']
        if (isinstance(phasecenter, str)
                and phasecenter.strip().upper() in ephem_sources):
            _ephemsrcname = phasecenter
        _imsize, _cell, _phasecenter = _handle_image_params(
            imsize, cell, phasecenter, sorted_vis,
            sorted_field, sorted_spw, sorted_antenna,
            sorted_scan, sorted_intent, sorted_timerange,
            _restfreq, pointingcolumn, _ephemsrcname
        )

        # Set up PySynthesisImager input parameters
        # - List all parameters that you need here
        # - Defaults will be assumed for unspecified parameters
        # - Nearly all parameters are identical to that in the task. Please look at the
        # list of parameters under __init__ using  " help ImagerParameters " )
        casalog.post('*** Creating paramList ***', origin=origin)
        paramList = ImagerParameters(
                # input file name
                msname=sorted_vis,  # 'sdimaging.ms',
                # data selection
                field=sorted_field,  # '',
                spw=sorted_spw,  # '0',
                timestr=sorted_timerange,
                antenna=sorted_baseline,
                scan=sorted_scan,
                state=sorted_intent,
                # image parameters
                imagename=output_path_prefix,  # 'try2',
                nchan=imnchan,  # 1024,
                start=imstart,  # '0',
                width=imwidth,  # '1',
                outframe=outframe,
                veltype=veltype,
                restfreq=_restfreq,
                phasecenter=_phasecenter,  # 'J2000 17:18:29 +59.31.23',
                imsize=_imsize,  # [75,75],
                cell=_cell,  # ['3arcmin', '3arcmin'],
                projection=projection,
                stokes=stokes,
                specmode=specmode,
                gridder='singledish',
                # single dish specific parameters
                gridfunction=gridfunction,
                convsupport=convsupport,
                truncate=gtruncate,
                gwidth=ggwidth,
                jwidth=gjwidth,
                pointingcolumntouse=pointingcolumn,
                convertfirst=convertfirst,
                minweight=minweight,
                clipminmax=clipminmax,
                # normalizer
                normtype='flatsky',
                pblimit=1e-16,  # TODO: explain why 1e-16 ?
                interpolation=interpolation,
                makesingledishnormalizer=True
            )

        # Construct the PySynthesisImager object, with all input parameters
        casalog.post('*** Creating imager object ***', origin=origin)
        imager = PySynthesisImager(params=paramList)

        # Initialize PySynthesisImager "modules"
        # required for single-dish imaging
        # - Pick only the modules you will need later on.
        # For example, to only make the PSF, there is no need for
        #  the deconvolver or iteration control modules.
        casalog.post('*** Initializing imagers ***', origin=origin)
        # This is where the underlying C++ synthesis imager is created
        imager.initializeImagers()
        casalog.post('*** Initializing normalizers ***', origin=origin)
        imager.initializeNormalizers()

        # Create Single-Dish images
        casalog.post('*** Creating single-dish images ***', origin=origin)
        imager.makeSdImage()
        casalog.post('*** Created single-dish images ***', origin=origin)

    finally:  # Close tools and rename Synthesis Imager's residual image
        casalog.post('*** Cleaning up tools ***', origin=origin)
        if imager is not None:
            imager.deleteTools()
        # Change image suffix from .residual to .image
        # residual_image_path = output_path_prefix + residual_suffix
        # if os.path.exists(residual_image_path):
        #     os.rename(residual_image_path, singledish_image_path)

    # Set single-dish image's beam size
    # TODO: re-define related functions in the new tool framework (sdms?)
    # ms_index = 0
    rep_ms = _get_param(0, infiles)
    rep_field = _get_param(0, field)
    rep_spw = _get_param(0, _spw)
    rep_antenna = _get_param(0, antenna)
    rep_scan = _get_param(0, scan)
    rep_intent = _get_param(0, intent)
    rep_timerange = _get_param(0, timerange)
    if len(rep_antenna) > 0:
        baseline = '{0}&&&'.format(rep_antenna)
    else:
        baseline = '*&&&'
    with open_ms(rep_ms) as ms:
        ms.msselect({'baseline': baseline})
        ndx = ms.msselectedindices()
        antenna_index = ndx['antenna1'][0]
    with sdutil.table_manager(os.path.join(rep_ms, 'ANTENNA')) as tb:
        antenna_name = tb.getcell('NAME', antenna_index)
        antenna_diameter = tb.getcell('DISH_DIAMETER', antenna_index)
    casalog.post("Setting single-dish image's beam")
    set_beam_size(
        rep_ms, singledish_image_path,
        rep_field, rep_spw, baseline, rep_scan, rep_intent, rep_timerange,
        _ephemsrcname, pointingcolumn, antenna_name, antenna_diameter,
        _restfreq, gridfunction, convsupport, truncate, gwidth, jwidth
    )

    # Set single-dish image's brightness unit (CAS-11503)
    if len(image_unit) == 0:
        image_unit = get_brightness_unit_from_ms(rep_ms)
    if len(image_unit) > 0:
        with open_ia(singledish_image_path) as ia:
            casalog.post(
                "Setting single-dish image's brightness unit to "
                f"'{image_unit}'")
            ia.setbrightnessunit(image_unit)

    # Update single-dish image's mask: mask low weight pixels
    weight_image_path = output_path_prefix + weight_suffix
    casalog.post("Creating weight image mask")
    do_weight_mask(singledish_image_path, weight_image_path, minweight)

    # Delete images systematically generated by the Synthesis Imager
    # which are either not required or currently useless
    # in the context of single-dish imaging
    # CAS-10891
    _remove_image(output_path_prefix + '.sumwt')

    # CAS-10893
    # TODO: remove the following line once the 'correct' SD
    # PSF image based on primary beam can be generated
    _remove_image(output_path_prefix + '.psf')
