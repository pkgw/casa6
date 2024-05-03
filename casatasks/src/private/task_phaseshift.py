from __future__ import absolute_import

from .mstools import write_history
from casatools import table, ms, mstransformer
from casatools import measures as me
from casatasks import casalog
from .parallel.parallel_data_helper import ParallelDataHelper


def phaseshift(
    vis=None, outputvis=None, keepmms=None, field=None,
    spw=None, scan=None, intent=None, array=None, 
    observation=None, datacolumn=None, phasecenter=None
):
    """
    Changes the phase center for either short or large
    offsets/angles w.r.t. the original
    """
    casalog.origin('phaseshift')

    if len(phasecenter) == 0:
        raise ValueError('phasecenter parameter must be specified')
    # Initiate the helper class
    pdh = ParallelDataHelper("phaseshift", locals())

    # Validate input and output parameters
    try:
        pdh.setupIO()
    except Exception as instance:
        casalog.post(str(instance), 'ERROR')
        raise RuntimeError(str(instance))

    # Input vis is an MMS
    if pdh.isMMSAndNotServer(vis) and keepmms:
        if not pdh.validateInputParams():
            raise RuntimeError('Unable to continue with MMS processing')

        pdh.setupCluster('phaseshift')

        # Execute the jobs
        try:
            pdh.go()
        except Exception as instance:
            casalog.post(str(instance), 'ERROR')
            raise RuntimeError(str(instance))
        return

    # Create local copies of tools (has to be here, otherwise
    # ParallelDataHelper has a problem digest the locals
    tblocal = table()
    mslocal = ms()
    mtlocal = mstransformer()

    # Actual task code starts here
    try:
        dirstr = phasecenter.split(' ')
        if not melocal.direction(dirstr[0], dirstr[1], dirstr[2]):
            raise ValueError("Illegal phacecenter specification " + phasecenter)
        try:
            # Gather all the parameters in a dictionary.
            config = {}

            config = pdh.setupParameters(
                inputms=vis, outputms=outputvis, field=field,
                spw=spw, array=array, scan=scan, intent=intent,
                observation=observation
            )

            # Check if CORRECTED column exists, when requested
            datacolumn = datacolumn.upper()
            if datacolumn == 'CORRECTED':
                tblocal.open(vis)
                if 'CORRECTED_DATA' not in tblocal.colnames():
                    casalog.post(
                        'Input CORRECTED_DATA does not exist. Will use DATA',
                        'WARN'
                    )
                    datacolumn = 'DATA'
                tblocal.close()

            casalog.post('Will use datacolumn = ' + datacolumn, 'DEBUG')
            config['datacolumn'] = datacolumn

            # Call MSTransform framework with tviphaseshift=True
            config['tviphaseshift'] = True
            config['reindex'] = False
            tviphaseshift_config = {'phasecenter': phasecenter}
            config['tviphaseshiftlib'] = tviphaseshift_config

            # Configure the tool
            casalog.post(str(config), 'DEBUG1')
            mtlocal.config(config)

            # Open the MS, select the data and configure the output
            mtlocal.open()

            # Run the tool
            casalog.post('Shift phase center')
            mtlocal.run()

        except Exception as instance:
            casalog.post(str(instance), 'ERROR')
            raise RuntimeError(str(instance))

        finally:
            mtlocal.done()

        # Write history to output MS, not the input ms.
        try:
            param_names = phaseshift.__code__.co_varnames[
                :phaseshift.__code__.co_argcount
            ]
            vars = locals()
            param_vals = [vars[p] for p in param_names]
            casalog.post('Updating the history in the output', 'DEBUG1')
            write_history(
                mslocal, outputvis, 'phaseshift', param_names,
                param_vals, casalog
            )
        except Exception as instance:
            casalog.post(
                "*** Error " + str(instance)
                + " updating HISTORY", 'WARN'
            )
            raise RuntimeError(str(instance))

        casalog.post('Updating the FIELD subtable of the output MeasurementSet with shifted'
                     ' phase centers', 'INFO')
        _update_field_subtable(outputvis, phasecenter)

    finally:
        if tblocal:
            tblocal.done()
            tblocal = None
        if mslocal:
            mslocal.done()
            mslocal = None
        if mtlocal:
            mtlocal.done()
            mtlocal = None


def _update_field_subtable(outputvis, phasecenter):
    """ Update MS/FIELD subtable with shifted center(s). """
    try:
        tblocal = table()
        # modify FIELD table
        tblocal.open(outputvis + '/FIELD', nomodify=False)
        pcol = tblocal.getcol('PHASE_DIR')

        if isinstance(phasecenter, str):
            thenewra_rad, thenewdec_rad = _convert_to_ra_dec_j2000(phasecenter)
            for row in range(0, tblocal.nrows()):
                pcol[0][0][row] = thenewra_rad
                pcol[1][0][row] = thenewdec_rad

        elif isinstance(phasecenter, dict):
            for field_id, field_center in phasecenter.items():
                thenewra_rad, thenewdec_rad = _convert_to_ra_dec_j2000(field_center)
                field_iidx = int(field_id)
                pcol[0][0][field_iidx] = thenewra_rad
                pcol[1][0][field_iidx] = thenewdec_rad

        tblocal.putcol('PHASE_DIR', pcol)

    except Exception as instance:
        casalog.post(
            "*** Error \'%s\' updating FIELD subtable" + str(instance),
            'WARN')
        raise RuntimeError(str(instance))
    finally:
        tblocal.done()


def _convert_to_ra_dec_j2000(phasecenter: str):
    """ Parse phase center string to obtain ra/dec (in rad) """
    dirstr = phasecenter.split(' ')
    try:
        melocal = me()
        thedir = melocal.direction(dirstr[0], dirstr[1], dirstr[2])
        if not thedir:
            raise RuntimeError(f"measures.direction() failed for phasecenter string:"
                               " {phasecenter}")
        if (dirstr[0] != 'J2000'):
            # Convert to J2000
            thedir = melocal.measure(thedir, 'J2000')
        thenewra_rad = thedir['m0']['value']
        thenewdec_rad = thedir['m1']['value']
    except Exception as instance:
        casalog.post(
            "*** Error " + str(instance)
            + " when interpreting parameter \'phasecenter\': ",
            'SEVERE'
        )
        raise RuntimeError(str(instance))
    finally:
        melocal.done()

    return thenewra_rad, thenewdec_rad
