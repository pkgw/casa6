from typing import Dict, Optional, Tuple, Union

import numpy as np
from .mstools import write_history
from casatools import table, ms, mstransformer
from casatools import measures as me
from casatasks import casalog
from .parallel.parallel_data_helper import ParallelDataHelper


def phaseshift(
    vis: str,
    outputvis: str,
    keepmms: bool,
    field: Optional[str],
    spw: Optional[str],
    scan: Optional[str],
    intent: Optional[str],
    array: Optional[str],
    observation: Optional[str],
    datacolumn: Optional[str],
    phasecenter: Union[str, dict],
):
    """
    Changes the phase center for either short or large
    offsets/angles w.r.t. the original
    """
    casalog.origin("phaseshift")

    if len(phasecenter) == 0:
        raise ValueError("phasecenter parameter must be specified")
    # Initiate the helper class
    pdh = ParallelDataHelper("phaseshift", locals())

    # Validate input and output parameters
    try:
        pdh.setupIO()
    except Exception as instance:
        casalog.post(str(instance), "ERROR")
        raise RuntimeError(str(instance))

    # Input vis is an MMS
    if pdh.isMMSAndNotServer(vis) and keepmms:
        if not pdh.validateInputParams():
            raise RuntimeError("Unable to continue with MMS processing")

        pdh.setupCluster("phaseshift")

        # Execute the jobs
        try:
            pdh.go()
        except Exception as instance:
            casalog.post(str(instance), "ERROR")
            raise RuntimeError(str(instance))
        return

    # Actual task code starts here
    # Gather all the parameters in a dictionary.
    config = {}

    config = pdh.setupParameters(
        inputms=vis,
        outputms=outputvis,
        field=field,
        spw=spw,
        array=array,
        scan=scan,
        intent=intent,
        observation=observation,
    )

    colnames = _get_col_names(vis)
    # Check if CORRECTED column exists, when requested
    datacolumn = datacolumn.upper()
    if datacolumn == "CORRECTED":
        if "CORRECTED_DATA" not in colnames:
            casalog.post(
                "Input data column CORRECTED_DATA does not exist. Will use DATA", "WARN"
            )
            datacolumn = "DATA"

    casalog.post(f"Will use datacolumn = {datacolumn}", "DEBUG")
    config["datacolumn"] = datacolumn

    # Call MSTransform framework with tviphaseshift=True
    config["tviphaseshift"] = True
    config["reindex"] = False
    tviphaseshift_config = {"phasecenter": phasecenter}
    config["tviphaseshiftlib"] = tviphaseshift_config

    # Configure the tool
    casalog.post(str(config), "DEBUG1")

    mtlocal = mstransformer()
    try:
        mtlocal.config(config)

        # Open the MS, select the data and configure the output
        mtlocal.open()

        # Run the tool
        casalog.post("Shift phase center")
        mtlocal.run()
    finally:
        mtlocal.done()

    # Write history to output MS, not the input ms.
    try:
        mslocal = ms()
        param_names = phaseshift.__code__.co_varnames[: phaseshift.__code__.co_argcount]
        vars = locals()
        param_vals = [vars[p] for p in param_names]
        casalog.post("Updating the history in the output", "DEBUG1")
        write_history(
            mslocal, outputvis, "phaseshift", param_names, param_vals, casalog
        )
    except Exception as instance:
        casalog.post(f"*** Error {instance} updating HISTORY", "WARN")
        raise RuntimeError(str(instance))
    finally:
        mslocal.done()

    casalog.post(
        "Updating the FIELD subtable of the output MeasurementSet with shifted"
        " phase centers",
        "INFO",
    )
    _update_field_subtable(outputvis, field, phasecenter)


def _get_col_names(vis: str) -> np.ndarray:
    tblocal = table()
    try:
        tblocal.open(vis)
        colnames = tblocal.colnames()
    finally:
        tblocal.done()
    return colnames


def _update_field_subtable(outputvis: str, field: str, phasecenter: Union[str, dict]):
    """Update MS/FIELD subtable with shifted center(s)."""

    try:
        tblocal = table()
        # modify FIELD table
        tblocal.open(outputvis + "/FIELD", nomodify=False)
        pcol = tblocal.getcol("PHASE_DIR")

        field_frames = _find_field_ref_frames(tblocal)
        if isinstance(phasecenter, str):
            if field:
                try:
                    field_id = int(field)
                except ValueError as _exc:
                    fnames = tblocal.getcol("NAME")
                    field_id = np.where(fnames == field)[0][0]
                thenewra_rad, thenewdec_rad = _convert_to_ref_frame(phasecenter, field_frames[field_id])
                pcol[0][0][field_id] = thenewra_rad
                pcol[1][0][field_id] = thenewdec_rad
            else:
                for row in range(0, tblocal.nrows()):
                    thenewra_rad, thenewdec_rad = _convert_to_ref_frame(phasecenter, field_frames[row])
                    pcol[0][0][row] = thenewra_rad
                    pcol[1][0][row] = thenewdec_rad

        elif isinstance(phasecenter, dict):
            for field_id, field_center in phasecenter.items():
                field_iidx = int(field_id)
                thenewra_rad, thenewdec_rad = _convert_to_ref_frame(field_center, field_frames[field_iidx])
                pcol[0][0][field_iidx] = thenewra_rad
                pcol[1][0][field_iidx] = thenewdec_rad

        tblocal.putcol("PHASE_DIR", pcol)

    except Exception as instance:
        casalog.post(f"*** Error '%s' updating FIELD subtable {instance}", "WARN")
        raise RuntimeError(str(instance))
    finally:
        tblocal.done()


def _find_field_ref_frames(tblocal: table) -> Dict[int, str]:
    """
    Given an open FIELD subtable, returns a dict of {field: reference_frame} for the PHASE_DIR
    column, where field is of type int, and reference_frame is of type str.

    This handles:
    - simple metadata where the reference frame is in the keywords of the PHASE_DIR column (same
      ref frame for all fields)
    - variable (per-field) reference frame metadata, usually in a PhaseDir_Ref additional column
    """

    dir_col = "PHASE_DIR"
    metainfo = tblocal.getcolkeyword(dir_col, "MEASINFO")
    nrows = tblocal.nrows()
    if "Ref" in metainfo:
        ref_frame = metainfo["Ref"]
        field_ref_frames = {field_id: ref_frame for field_id in range(0, nrows)}
    elif "Ref" not in metainfo and ("VarRefCol" in metainfo and "TabRefTypes" in metainfo and
                                    "TabRefCodes" in metainfo):
        col = metainfo["VarRefCol"]   # presumably PhaseDir_Ref
        ref_frame_codes = tblocal.getcol(col)
        ref_codes_to_frame_idx = {code: idx for idx, code in enumerate(metainfo["TabRefCodes"])}
        ref_frame_names = [metainfo["TabRefTypes"][ref_codes_to_frame_idx[code]] for code in
                           ref_frame_codes]
        field_ref_frames = {field_id: ref_frame_names[field_id] for field_id in np.arange(0, nrows)}
    else:
        raise RuntimeError("Error when retrieving reference frames from the metadata of column "
                           f"{dir_col}. The field 'Ref' is not present but could not find the "
                           "fields 'VarRefCol', 'TabRefTypes', and 'TabRefCodes'")

    return field_ref_frames


def _convert_to_ref_frame(phasecenter: str, output_ref_frame: str) -> Tuple[float, float]:
    """
    Converts one phasecenter (presumably given as input to this task) to another reference
    frame (presumably the frame used in the (output) MS for the relevant field).

    When applying phaseshift, the output MS has the same reference frames as the input MS, as
    propagated by mstransform.

    Returns the v0,v1 (RA,Dec) values in units of radians and using the requested frame, ready to be
    written to a FIELD subtable.
    """
    def parse_phasecenter(center: str) -> Tuple[str, str, str]:
        """
        Parse phase center string to obtain ra/dec (in rad). Splits the:
        - (optional) frame,
        - v0 (typically RA),
        - v1 (typically Dec)
        in a phasecenter string.
        """
        dirstr = center.split(" ")
        if 3 == len(dirstr):
            dir_frame, dir_v0, dir_v1 = dirstr[0], dirstr[1], dirstr[2]
        elif 2 == len(dirstr):
            dir_frame = ""  # J2000 will be default in me.direction, etc.
            dir_v0, dir_v1 = dirstr[0], dirstr[1]
        else:
            raise AttributeError(f"Wrong phasecenter string: '{center}'. It must have 3 or 2 "
                                 "items separated by space(s): 'reference_frame RA Dec' or 'RA Dec'")

        return dir_frame, dir_v0, dir_v1

    try:
        melocal = me()
        dir_frame, dir_v0, dir_v1 = parse_phasecenter(phasecenter)

        # Note: even if the input frame is the same as output_ref_frame, we need to ensure units of
        # radians for the FIELD subtable
        if dir_frame:
            thedir = melocal.direction(dir_frame, dir_v0, dir_v1)
        else:
            thedir = melocal.direction(v0=dir_v0, v1=dir_v1)
        if not thedir:
            raise RuntimeError(
                f"measures.direction() failed for phasecenter string: {phasecenter}"
            )

        if dir_frame != output_ref_frame:
            thedir = melocal.measure(thedir, output_ref_frame)
        thenewra_rad = thedir["m0"]["value"]
        thenewdec_rad = thedir["m1"]["value"]
    except Exception as instance:
        casalog.post(
            f"*** Error {instance} when interpreting parameter 'phasecenter': ",
            "SEVERE",
        )
        raise RuntimeError(str(instance))
    finally:
        melocal.done()

    return thenewra_rad, thenewdec_rad
