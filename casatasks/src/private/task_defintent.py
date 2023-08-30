from __future__ import absolute_import
import os
import numpy as np
import shutil
import casatools

def defintent(vis='', intent='', mode='append',
              outputvis='', scan='', field='',
              obsid=''):
    """
    Description:
    Allows users to manually set the intents for a selection of scans, fields, or obsids.
    
    Keyword arguments:
    mode: 'set' or 'append'.
        Set allows the user to fully define a new intent.
        Append allows the user to add to the current intent for the current selection.
        Accepts string values, no change if left undefined.
    intent: new intent to add.
        User provides the new intent to be set.
        Accepts string values, no change if left undefined.
    scan: Select the scans to modify.
        Defaults to all scans selected.
    field: Select the fields to modify
        Defaults to all fields selected.
    obsid: Select the obsids to modify
        Defaults to all obsids selected
        
    Return: none
    """

    tb = casatools.table()
    ms = casatools.ms()
    #changeList = []
    
    # If no intent has been provided exit the task and print
    if vis == '':
        print('You must specify a MS')
        return
        
    if intent == '':
        print('you must specify an Intent')
        return
        
    # if there is an outputvis make a copy
    if outputvis != '' and outputvis != vis:
        if os.path.exists(outputvis):
            print("outputvis already exists! Exiting task...")
            return
        shutil.copytree(vis, outputvis)
    # if the outputvis is the same as the vis. edit the vis table and don't make a copy
    elif outputvis == vis:
        outputvis = vis
        print("outputvis and vis are the same. Editing provided vis...")
    else:
        print("No outputvis has been specified, please enter an outputvis name")
        return
    
    # Table tool query?
    # ----- TABLE SELECTION -----
    
    # Get field names
    tb.open(vis+'/FIELD')
    fieldnames = tb.getcol('NAME')
    tb.close()
    
    # Get field ids for names
    fieldSplit = field.split(',')
    fieldSplit = [x.strip() for x in fieldSplit]
    
    # Link field names to ID
    nameDict = {}
    for i in range(len(fieldnames)):
        if fieldnames[i] not in nameDict:
            nameDict[fieldnames[i]] = str(i)
        else:
            nameDict[fieldnames[i]] += ',' + str(i)
            
    # Replace names in selection with IDs
    for i in range(len(fieldSplit)):
        if fieldSplit[i] in fieldnames:
            field = field.replace(fieldSplit[i], nameDict[fieldSplit[i]])
            
    # Get exitsing intents
    tb.open(vis+'/STATE')
    intentcol = tb.getcol('OBS_MODE')
    tb.close()
    
    # Allowed Intents? / Intents reformatting?
    
    # If selected field is found
    foundField = False
    if (type(scan) != list):
        scan = str(scan)
    
    #selectedRows = set()
    selectedRows = []
    selectedIntents = dict()
    
    # NEW get query using ms tool selection
    ms.open(vis)
    ms.msselect({'field':field, 'scan':scan, 'observation':obsid}, onlyparse=True)
    selectedIndex = ms.msselectedindices()
    ms.close()
    
    tb.open(vis)
    fieldIds = tb.getcol('FIELD_ID')
    scanNum = tb.getcol('SCAN_NUMBER')
    stateIds = tb.getcol('STATE_ID')
    obsIds = tb.getcol('OBSERVATION_ID')
    
    '''# Dict to write to the outfile if it exists
    outfileDict = {}
    outfileDict["origin_state_ids"] = stateIds
    paramDict = {'vis':vis, 'intent':intent, 'mode':mode, 'outfile':outfile,
              'originfile':originfile, 'scan':scan, 'field':field,
              'obsid':obsid}
    outfileDict["task_parameters"] = paramDict
    outfileDict["execution_time"] = date.today().strftime("%B %d, %Y")'''
    
    
    # mstool query version
    toJoin = []
    if len(selectedIndex['field']) > 0:
        toJoin.append(f"FIELD_ID in {list(selectedIndex['field'])}")
    if len(selectedIndex['scan']) > 0:
        toJoin.append(f"SCAN_NUMBER in {list(selectedIndex['scan'])}")
    if len(selectedIndex['observationid']) > 0:
        toJoin.append(f"OBSERVATION_ID in {list(selectedIndex['observationid'])}")
    # join into query string
    taskQuery = " && ".join(toJoin)
    
    selectedData = tb.query(taskQuery)
    #selectedRows = set(selectedData.rownumbers())
    selectedRows = selectedData.rownumbers()
    
    selectedStateIds = selectedData.getcol('STATE_ID')
    for i in range(len(selectedRows)):
        selectedIntents[selectedStateIds[i]] = selectedStateIds[i]
        
        #tmpString =  str(selectedRows[i]) + ':' + str(selectedStateIds[i])
        #changeList.append(tmpString)
        
    tb.close()
                
    print("Number of matching rows found: ", len(selectedRows))
    
    # for Set if intent not in state table
    # then add a new row to the state table and change index (STATE_ID) in main table
    if mode.lower() == 'set':
        # Keep track of the new value to set the state_id to
        newState = -1
        # Adding to intents col
        statetb = outputvis+'/STATE'
        tb.open(statetb, nomodify=False)
        intents = tb.getcol('OBS_MODE')
        # Check if the intent already exists
        if intent in intents:
            print("Intent already exists")
            newState = np.where(intents == intent)
        # If it doesn't add a row with the new intent
        else:
            # If there are no intents to begin with add blank zero to prevent listobs segfault (?)
            numIntents = len(list(tb.getcol('OBS_MODE')))
            if numIntents == 0:
                tb.addrows(2)
            else:
                tb.addrows(1)
            intents = list(tb.getcol('OBS_MODE'))
            intents[-1] = intent
            if numIntents == 0:
                print("No intents present. filling unselected with UNSPECIFIED(DEFINTENT)")
                intents[0] = 'UNSPECIFIED(DEFINTENT)'
            intents = np.asarray(intents)
            tb.putcol('OBS_MODE', intents)
            newState = len(intents) - 1
            tb.close()
            
            # For all selected rows replace with new state_id
            tb.open(outputvis, nomodify=False)
            stateCol = tb.getcol('STATE_ID')
            
            for row in selectedRows:
                stateCol[row] = newState
            tb.putcol('STATE_ID', stateCol)
            tb.close()
            
        tb.close()
    
    # For Append mode
    elif mode.lower() == 'append':
        statetb = outputvis+'/STATE'
        # Find our selected intents
        for i in selectedIntents:
            newState = -1
            tb.open(statetb, nomodify=False)
            intents = tb.getcol('OBS_MODE')
            # Add a row with old intent + new, if it was UNSPECIFIED just do the new intent (?)
            if intents[i] == 'UNSPECIFIED(DEFINTENT)':
                newIntent = intent
            else:
                newIntent = intents[i] + ',' + intent
            # Check if thie intent already exists
            if newIntent in intents:
                print("Intent already exists")
                newState = intents.index(newIntent)
            else:
                tb.addrows(1)
                intents = list(tb.getcol('OBS_MODE'))
                intents[-1] = newIntent
                intents = np.asarray(intents)
                tb.putcol('OBS_MODE', intents)
                newState = len(intents) - 1
                selectedIntents[i] = newState
                tb.close()
            tb.close()
        
        # For all selected rows replace with new ID
        tb.open(outputvis, nomodify=False)
        stateCol = tb.getcol('STATE_ID')
        for row in selectedRows:
            if stateCol[row] in selectedIntents:
                stateCol[row] = selectedIntents[stateCol[row]]
        tb.putcol('STATE_ID', stateCol)
        tb.close()

    return
