import os
import sys
import shutil

from casatasks import casalog
from casatools import ms as mstool
from casatools import table as tbtool
from .mslisthelper import check_mslist

def testconcat(vislist,testconcatvis,freqtol,dirtol,copypointing):
    """
    The list of data sets given in the vis argument are concatenated into an output
    data set in concatvis without copying the bulk data of the main table.
    This is useful for obtaining the information in the merged subtables without
    actually performing a time-consuming concatenation disk.

    Keyword arguments:
    vis -- Name of input visibility files for which the subtables are to be combined
    default: none; example: vis = 'mydata.ms',
    example: vis=['src2.ms','ngc5921.ms','ngc315.ms']
    
    testconcatvis -- Name of MS that will contain the concatenated subtables
    default: none; example: testconcatvis='test.ms'
    
    freqtol -- Frequency shift tolerance for considering data to be in the same
    spwid.  The number of channels must also be the same.
    default: ''  do not combine unless frequencies are equal
    example: freqtol='10MHz' will not combine spwid unless they are
    within 10 MHz.
        Note: This option is useful to conbine spectral windows with very slight
    frequency differences caused by Doppler tracking, for example.
    
    dirtol -- Direction shift tolerance for considering data as the same field
    default: '' means always combine.
    example: dirtol='1.arcsec' will not combine data for a field unless
    their phase center differ by less than 1 arcsec.  If the field names
    are different in the input data sets, the name in the output data
    set will be the first relevant data set in the list.

    copypointing -- copy all rows of the pointing table
    default: True
    """

    ###
    try:
        casalog.origin('testconcat')
        t = tbtool()
        m = mstool()
        #break the reference between vis and vislist as we modify vis
        if(type(vislist)==str):
            vis=[vislist]
        else:
            vis=list(vislist)

            # test the consistency of the setup of the different MSs
            casalog.post('Checking MS setup consistency ...', 'INFO')
            try:
                    mydiff = check_mslist(vis, ignore_tables=['SORTED_TABLE', 'ASDM*'], testcontent=False) 
            except Exception as instance:
                    raise RuntimeError("*** Error \'%s\' while checking MS setup consistency" % (instance))

            if mydiff != {}:
                    casalog.post('The setup of the input MSs is not fully consistent. The concatenation may fail', 'WARN')
                    casalog.post('and/or the affected columns may contain partially only default data.', 'WARN')
                    casalog.post(str(mydiff), 'WARN')

        if((type(testconcatvis)!=str) or (len(testconcatvis.split()) < 1)):
            raise Exception ('parameter testconcatvis is invalid')
        if(vis.count(testconcatvis) > 0):
            vis.remove(testconcatvis)

        if(os.path.exists(testconcatvis)):
            raise Exception ('Visibility data set '+testconcatvis+' exists. Will not overwrite.')
        else:
            if(len(vis) >0): 
                casalog.post('copying structure of '+vis[0]+' to '+testconcatvis , 'INFO')
                # Copy procedure which does not copy the bulk data of the main table
                t.open(vis[0])
                tmptb = t.copy(newtablename=testconcatvis, deep=True, valuecopy=True, norows=True) # copy only table structure
                tmptb.close()
                t.close()
                # copy content of subtables
                thesubtables = [f.name for f in os.scandir(vis[0]) if f.is_dir()]

                for subt in thesubtables:
                    if not (subt[0]=='.'):
                        t.open(vis[0]+'/'+subt)
                        no_rows = False
                        if (subt=='POINTING' and not copypointing):
                            casalog.post('*** copypointing==False: resulting MS will have empty POINTING table', 'INFO')
                            no_rows = True
                        tmptb = t.copy(testconcatvis+'/'+subt, deep=False, valuecopy=True, norows=no_rows)
                        tmptb.close()
                        t.close()
                vis.remove(vis[0])
        # determine handling switch value
        handlingswitch = 1
        if not copypointing:
            handlingswitch = 3
    
        m.open(testconcatvis,False) # nomodify=False to enable writing
    
        for elvis in vis : 
            casalog.post('concatenating subtables from '+elvis+' into '+testconcatvis , 'INFO')

            m.concatenate(msfile=elvis,freqtol=freqtol,dirtol=dirtol,handling=handlingswitch)

            m.writehistory(message='taskname=testconcat',origin='testconcat')
            m.writehistory(message='vis         = "'+str(testconcatvis)+'"',origin='testconcat')
            m.writehistory(message='concatvis   = "'+str(elvis)+'"',origin='testconcat')
            m.writehistory(message='freqtol     = "'+str(freqtol)+'"',origin='testconcat')
            m.writehistory(message='dirtol      = "'+str(dirtol)+'"',origin='testconcat')
            m.writehistory(message='copypointing = "'+str(copypointing)+'"',origin='concat')

        m.close()

    except Exception as instance:
        print('*** Error ***',instance)
        raise Exception (instance)

