
import os
import math
import numpy as np
import fnmatch

import subprocess
from collections import OrderedDict as odict

###some helper tools
from casatasks import casalog
from casatools import table, msmetadata
from casatools import quanta as qatool
from casatools import ms as mstool

ms = mstool()
tb = table()
msmd = msmetadata()


def check_mslist(vis, ignore_tables=['SORTED_TABLE'], testcontent=True):
    """
    Check the consistency of the setup of the MSs in the list "vis"
    Returns a dictionary describing the inconsistencies w.r.t. to the first MS in the list:
       {'<vis 1>':
              '<tablename1>': {'present_a': True/False,
                  'present_b': True/False,
                  'missingcol_a':[<column1>, <column2>,...],
                  'missingcol_b':[<column1>, <column2>,...]},
              '<tablename2>: {'present_a': True/False,
                  'present_b': True/False,
                  'missingcol_a':[<column1>, <column2>,...],
                  'missingcol_b':[<column1>, <column2>,...]}

         '<vis 2>':
                ...
        }

    where <vis n> is the name of the nth MS in the input list of MSs. An entry for a given MS
    is only present if there are differences between that MS and the first MS in the list.

    If there are no differences in the setup of the MSs, the returned dictionary is empty.

    If there are differences for an MS, there is a dictionary item for each table which has a different
    setup and the value is a dictionary with two lists of the names of the columns which are
    missing in the table of MS A (the first one in the list) and table of MS B (the one compared to).
    "Missing" is to be understood as "present in the other table but not in this one".
    Furthermore, the dictionary contains the items "present_a" and "present_b" which are True
    if the given table is present at all in MS A and MS B respectively.

    The optional parameter "ignore_tables" defines the list of subtables of Main which
    are to be ignored in the comparison. Default: ['SORTED_TABLE']
    Table names can be provided using wildcards like "*" and "?", e.g. 'ASDM_*'.

    If the optional parameter testcontent==True, then for a column which is absent in one table
    it is tested in the other table whether the column actually contains data,
    i.e. cell 0 can be read. If not, the absence of the column is ignored.

    Independently from the value of testcontent, all optional Main table columns are
    tested as to whether they are present and if so whether they contain data.
    A warning is raised if they don't contain data.

    """

    rval = {}

    if type(vis) != list:
        if type(vis) == str:
            vis = [vis]
        else:
            raise ValueError('vis parameter needs to be a list of strings.')

    if type(ignore_tables) != list:
        if type(ignore_tables) == str:
            ignore_tables = [ignore_tables]
        else:
            raise ValueError('ignore_tables parameter needs to be a list of strings.')


    if len(vis) == 1:
        try:
            ms.open(vis[0])
            ms.close()
        except:
            raise ValueError(vis[0]+' does not exist or is not a MeasurementSet.')

        return rval

    haspointing = np.zeros(len(vis)) # track the presence of pointing tables
    viscount = 0

    # Gather information from first MS in list

    tb.open(vis[0])
    descr_a = tb.getdesc()
    tb.close()

    descr_a['_name_'] = vis[0]

    descr_a_kw = descr_a['_keywords_']
    if not 'MS_VERSION' in descr_a_kw:
        raise ValueError(vis[0]+' is not a MeasurementSet.')

    # Eliminate the tables to be ignored
    tbdel = []
    for mytablepattern in ignore_tables:
        for mytable in descr_a_kw:
            if fnmatch.fnmatch(mytable, mytablepattern):
                tbdel.append(mytable)
    for mytable in tbdel:
        del descr_a_kw[mytable]

    # Extract subtable details
    subtbpaths_a = []
    subtbnames_a = []
    subtbdescs_a = []

    for mysubtb in descr_a_kw:
        if type(descr_a_kw[mysubtb]) == str:
            subtbpath = descr_a_kw[mysubtb].split(' ')
            if subtbpath[0] == 'Table:':
                subtbpaths_a.append(subtbpath[1])
                myname = subtbpath[1].split('/')[-1]
                subtbnames_a.append(myname)
                tb.open(subtbpath[1])
                mydesc = tb.getdesc()
                if myname == 'POINTING':
                    haspointing[0] = 1
                    casalog.post('Checking for unpopulated POINTING table in first MS ...', 'INFO')
                    try:
                        tb.getcell('TIME',0)
                    except:
                        haspointing[0] = 0
                tb.close()
                mydesc['_name_'] = subtbpath[1]
                subtbdescs_a.append(mydesc)

    casalog.post('Checking for unpopulated optional Main Table columns in first MS ...', 'INFO')
    opt_main_populated(descr_a) # ... in first MS

    # Loop over other MSs and check against first

    for myvis in vis[1:]:
        if myvis==vis[0]:
            raise ValueError(myvis+' is contained in the list more than once.')

        viscount += 1

        tb.open(myvis)
        descr_b = tb.getdesc()
        tb.close()

        descr_b['_name_'] = myvis

        descr_b_kw = descr_b['_keywords_']
        if not 'MS_VERSION' in descr_b_kw:
            raise ValueError(myvis+' is not a MeasurementSet.')

        # Eliminate the tables to be ignored
        tbdel = []
        for mytablepattern in ignore_tables:
            for mytable in descr_b_kw:
                if fnmatch.fnmatch(mytable, mytablepattern):
                    tbdel.append(mytable)
        for mytable in tbdel:
            del descr_b_kw[mytable]

        # Extract subtable details
        subtbpaths_b = []
        subtbnames_b = []
        subtbdescs_b = []

        for mysubtb in descr_b_kw:
            if type(descr_b_kw[mysubtb]) == str:
                subtbpath = descr_b_kw[mysubtb].split(' ')
                if subtbpath[0] == 'Table:':
                    subtbpaths_b.append(subtbpath[1])
                    myname = subtbpath[1].split('/')[-1]
                    subtbnames_b.append(myname)
                    tb.open(subtbpath[1])
                    mydesc = tb.getdesc()
                    if myname == 'POINTING':
                        haspointing[viscount] = 1
                        casalog.post('Checking for unpopulated POINTING table ...', 'INFO')
                        try:
                            tb.getcell('TIME',0)
                        except:
                            haspointing[viscount] = 0
                    tb.close()
                    mydesc['_name_'] = subtbpath[1]
                    subtbdescs_b.append(mydesc)

        # Comparison
        compresult = {}

        # Main table
        cmpres = comptbdescr(descr_a, descr_b, testcontent=testcontent)
        if cmpres != {}:
            compresult['Main'] = cmpres

        casalog.post('Checking for unpopulated optional Main Table columns ...', 'INFO')
        opt_main_populated(descr_b)

        # Subtables
        for i in range(len(subtbnames_a)): # loop over tables in first MS
            if not subtbnames_a[i] in subtbnames_b:
                compresult[subtbnames_a[i]] = {'present_a': True, 'present_b': False}
            else: # table is present in both MSs
                cmpres = comptbdescr(subtbdescs_a[i], subtbdescs_b[ subtbnames_b.index(subtbnames_a[i]) ],
                                     testcontent=testcontent)
                if cmpres != {}:
                    compresult[subtbnames_a[i]] = cmpres

        for i in range(len(subtbnames_b)): # loop over tables in second MS
            if not subtbnames_b[i] in subtbnames_a:
                compresult[subtbnames_b[i]] = {'present_a': False, 'present_b': True}
            # else clause not needed since already covered in previous loop

        if compresult != {}:
            rval[myvis] = compresult

    # evaluate haspointing array
    if (1 in haspointing) and (False in ( haspointing == 1 )):
        casalog.post('Some but not all of the input MSs are lacking a populated POINTING table:', 'WARN')
        for i in range(len(haspointing)):
            if haspointing[i] == 0:
                casalog.post('   '+str(i)+': '+vis[i], 'WARN')
        casalog.post('The joint dataset will not have a valid POINTING table.', 'WARN')

    return rval


def comptbdescr(descr_a, descr_b, ignorecol=[], testcontent=True):
    """Utility function for check_mslist
       - compares two table descriptions descr_a and descr_b
       - the absence of the columns listed in ignorecol is ignored
       - if testcontent==True, then for a column which is absent in one table
         it is tested in the other table whether the column actually contains data,
         i.e. cell 0 can be read. If not, the absence of the column is ignored.
         For this to work, the table path has to be added to the table description
         as item "_name_".
    """
    rval = {}
    mscol_a = []
    mscol_b = []
    for myentry in descr_a:
        if myentry[0]!='_' and not myentry in ignorecol: # only inspect relevant columns
            if not myentry in descr_b:
                if testcontent:
                    tb.open(descr_a['_name_'])
                    try:
                        tb.getcell(myentry,0)
                    except:
                        tb.close()
                        casalog.post('Column '+myentry+' in table '+descr_a['_name_']+' has no data.','INFO')
                        continue # i.e. ignore this column because it has no data
                    tb.close()
                mscol_b.append(myentry)
    for myentry in descr_b:
        if myentry[0]!='_' and not myentry in ignorecol: # only inspect relevant columns
            if not myentry in descr_a:
                if testcontent:
                    tb.open(descr_b['_name_'])
                    try:
                        tb.getcell(myentry,0)
                    except:
                        tb.close()
                        casalog.post('Column '+myentry+' in table '+descr_b['_name_']+' has no data.','INFO')
                        continue # i.e. ignore this column because it has no data
                    tb.close()
                mscol_a.append(myentry)
    if mscol_a!=[] or mscol_b!=[]:
        rval = {'present_a': True, 'present_b': True,
                'missingcol_a': mscol_a, 'missingcol_b': mscol_b}

    return rval

def sort_mslist(vis, visweightscale=None):
    """
    Returns two or three items:
      1) list of MSs sorted by the earliest entry in the Main table TIME column.
      2) list of sorted MS start times
      3) if visweightscale!=None and contains a list of corresponding numbers,
         they are sorted as well and returned as third return value.
         If visweightscale==[], a list filled with values of 1 is returned.

    vis - list of MS names
    visweightscale - list of numbers (e.g. the weight scaling factors in concat)
               default: None (no value provided)

    """
    if type(vis) != list:
        if type(vis)==str:
            vis = [vis]
        else:
            raise ValueError('Parameter vis should be a list of strings.')

    doweightscale = True
    if type(visweightscale)!=list:
        if visweightscale!=None:
            try:
                visweightscale = [float(visweightscale)]
            except:
                raise ValueError('Parameter visweightscale should be a list of numerical values or None.')
        else:
            doweightscale = False
    elif visweightscale==[]:
        visweightscale = list(np.ones(len(vis)))
    elif len(visweightscale) != len(vis):
        raise ValueError('Parameter visweightscale should have same length as vis.')


    sortedvis = []
    sortedvisweightscale = []
    sortedtimes = []
    namestuples = []
    for name in vis:
        tb.open(name)
        times = tb.getcol('TIME')
        tb.close()
        times.sort()
        if doweightscale:
            namestuples.append( (times[0], name, visweightscale[vis.index(name)]) )
        else:
            namestuples.append( (times[0], name, 0) )

    sorted_namestuples = sorted(namestuples, key=lambda msname: msname[0])

    for i in range(0,len(vis)):
        sortedvis.append(sorted_namestuples[i][1])
        sortedtimes.append(sorted_namestuples[i][0])
        if doweightscale:
            sortedvisweightscale.append(sorted_namestuples[i][2])

    if doweightscale:
        return sortedvis, sortedtimes, sortedvisweightscale
    else:
        return sortedvis, sortedtimes


def report_sort_result(sorted_vis, sorted_times, sorted_idx, mycasalog=None, priority='INFO'):
    """Report result of MS sort.

    Args:
        sorted_vis (list): sorted list of MS
        sorted_times (list): sorted list of observation start time
        sorted_idx (list): list of indices of original order of MS list
        mycasalog (logsink, optional): logsink instance for logging. Defaults to None.
        priority (str, optional): priority for logging. Defaults to 'WARN'.
    """
    if len(sorted_vis) <= 1:
        # trivial result. do nothing.
        return

    if mycasalog is None:
        local_casalog = casalog
    else:
        local_casalog = mycasalog
    qa = qatool()
    header = 'Order {:>24s} {:>20s} Original_Order'.format('MS_Name', 'Start_Time')
    local_casalog.post('Summary of the MS internal sort:', priority=priority)
    local_casalog.post(header, priority=priority)
    local_casalog.post('-' * len(header), priority=priority)
    for isort, (iorig, v, t) in enumerate(zip(sorted_idx, sorted_vis, sorted_times)):
        local_casalog.post(
            '{:>3d} {:>26s} {:>20s} {:>3d}'.format(
                isort,
                os.path.basename(v.rstrip('/')),
                qa.time(qa.quantity(t, 's'), form=['ymd', 'hms'])[0],
                iorig
            ),
            priority=priority
        )


def opt_main_populated(descr, ignorecol=[]):
    """Utilty function for check_mslist
       Check the optional Main Table data columns and raise warnings
       if they exist but don't contain data.

       descr - table description of the main table

       The absence of columns listed in ignorecol is ignored.

       Returns True if no warnings were raised.
    """

    rval = True

    opt_main_cols = ['DATA', 'FLOAT_DATA', 'LAG_DATA', 'SIGMA_SPECTRUM', 'WEIGHT_SPECTRUM']

    tbname = descr['_name_']

    for myentry in opt_main_cols:
        if myentry in descr and not myentry in ignorecol: # only inspect relevant columns
            tb.open(tbname)
            try:
                tb.getcell(myentry,0)
            except:
                tb.close()
                rval = False
                casalog.post('Column '+myentry+' in table '+tbname+' has no data. Accessing it will cause errors.','WARN')
                continue
            tb.close()

    return rval



