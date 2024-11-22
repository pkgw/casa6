import json

from casatools import table
from casatasks import casalog

from .correct_ant_posns_evla import correct_ant_posns_evla as _correct_ant_posns_evla

_tb = table( )

def correct_ant_posns(vis_name, print_offsets=False, time_limit=0):
    """
    Given an input visibility MS name (vis_name), find the antenna
    position offsets that should be applied.  This application should
    be via the gencal task, using caltype='antpos'.

    If the print_offsets parameter is True, will print out each of
    the found offsets (or indicate that none were found), otherwise
    runs silently.

    A list is returned where the first element is the returned error
    code, the second element is a string of the antennas, and the
    third element is a list of antenna Bx,By,Bz offsets. An example
    return list might look like:
    [ 0, 'ea01,ea19', [0.0184, -0.0065, 0.005, 0.0365, -0.0435, 0.0543] ]

    The second and third elements of the list returned are in the format
    expected by the calibrater tool method cb.specifycal() for the
    parameters antenna and parameter, respectively.

    Usage examples:

       CASA <1>: antenna_offsets = correct_ant_posns('test.ms')
       CASA <2>: if (antenna_offsets[0] == 0):
       CASA <3>:     gencal(vis='test.ms', caltable='cal.G', \
                     caltype='antpos', antenna=antenna_offsets[1], \
                     parameter=antenna_offsets[2])

    For specific details for the EVLA see correct_ant_posns_evla.
    """
    _tb.open(vis_name+'/OBSERVATION')
    # specific code for different telescopes
    tel_name = _tb.getcol('TELESCOPE_NAME')
    _tb.close()
    if tel_name == 'EVLA' or tel_name == 'VLA':
        return _correct_ant_posns_evla(vis_name, print_offsets, time_limit)
    else:
        msg = 'Currently only work for EVLA observations'
        if (print_offsets):
            print(msg)
        else:
            # send to casalogger
            casalog.post(msg, "WARN")
        return [1, '', []]

def correct_ant_posns_alma_json(vis_name, json_file):
    """
    This function computes the differences of the ALMA antenna
    positions given the corrected positions in a JSON file (presumably
    created by the getantposalma task).

    For each antenna found in the JSON file which is also present in the
    ANTENNA subtable, the difference between the corrected position in
    the JSON file and the nominal position in the MS is computed.

    If the MS does not belong to an ALMA observation an exception is
    thrown.

    Likewise, if the JSON file has not been created by the getantposalma an
    exception is thrown.

    It returns a tuple with the antenna names as a string with the antenna
    names separated by commma and the differential possitions as a long
    vector with a list of antenna Bx,By,Bz offsets. An example return
    might look like:
    ['DA42,DA43', [-0.001, 0.001, 0.001, 0.002, -0.002, 0.002] ]

    These two elements of the list returned are in the format
    expected by the calibrater tool method cb.specifycal() for the
    parameters antenna and parameter, respectively.
    """

    # get the telescope name
    _tb.open(vis_name+'/OBSERVATION')
    tel_name = _tb.getcol('TELESCOPE_NAME')
    _tb.close()
    with open(json_file, "r") as f:
        corrected_antenna_abspos_map = json.load(f)

    # throw if the MS wasn't observed with ALMA
    if tel_name[0] != 'ALMA' :
        raise ValueError('Antenna positions are from ALMA but MS is from '+tel_name[0]+' telescope')

    # throw if the JSON file wasn't generated using the ALMA web service
    if(corrected_antenna_abspos_map['metadata']['product_code'] != 'antposalma') :
        raise ValueError('JSON file with antenna positions is not from ALMA. This is currently not supported')

    # compute the differences in antenna positions between the corrected ones in the
    # JSON file and the nominal ones in the MS.
    _tb.open(vis_name+'/ANTENNA')
    nominal_antenna_names = list(_tb.getcol('NAME'))
    nominal_abspos = _tb.getcol('POSITION')
    antennas_to_correct = ""
    differential_pos = []
    for corrected_antenna_name, corrected_abspos in corrected_antenna_abspos_map['data'].items() :
        if corrected_antenna_name in nominal_antenna_names :
            tb_idx = nominal_antenna_names.index(corrected_antenna_name)
            if len(antennas_to_correct) != 0 :
                antennas_to_correct += ","
            antennas_to_correct += corrected_antenna_name
            differential_pos.append(corrected_abspos[0] - nominal_abspos[0][tb_idx])
            differential_pos.append(corrected_abspos[1] - nominal_abspos[1][tb_idx])
            differential_pos.append(corrected_abspos[2] - nominal_abspos[2][tb_idx])
        else:
            raise ValueError("Antenna {} could not be found in input MS".format(corrected_antenna_name))
    if (len(differential_pos) == 0) :
        raise ValueError("The list of antenna positions in JSON file is empty")


    return antennas_to_correct, differential_pos

