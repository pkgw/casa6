from __future__ import absolute_import
import os
import numpy as np

import casatools

def defintent(vis='', intent='', mode='',
              scan='', field='', obsid='',
              revertList=[]):
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
        
    Return: Returns a list detailing which STATE_IDs have been changed. This can be used to revert changes
    """

    tb = casatools.table()
    ms = casatools.ms()
    changeList = []
    
    # If no intent has been provided exit the task and print
    if vis == '':
        print('You must specify a MS')
        return
        
    if intent == '':
        print('you must specify an Intent')
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
        
        tmpString =  str(selectedRows[i]) + ':' + str(selectedStateIds[i])
        changeList.append(tmpString)
        
    tb.close()
                
    print("Number of matching rows found: ", len(selectedRows))
    print(mode.lower())
    
    # For Revert mode
    if mode.lower() == 'revert':
        # If there is no provided revertList break out
        if revertList == []:
            print("No revert list provided")
            
        # Revert list structure is [(<ROW>,<OLD STATE_ID>)...]
        # Iterate over all the changed rows and set the state id to the old one
        tb.open(vis, nomodify=False)
        stateCol = tb.getcol('STATE_ID')
        #rowsToRemove = set()
        for item in revertList:
            stateCol[int(item.split(':')[0])] = int(item.split(':')[1])
            #rowsToRemove.add(int(item.split(':')[1]))
        # Set the state col back after reverting
        tb.putcol('STATE_ID', stateCol)
        tb.close()
        
        # Remove the row from the STATE table
        tb.open(vis+'/STATE', nomodify=False)
        modes = tb.getcol('OBS_MODE')
        for i in range(len(modes)):
            if modes[i] == intent:
                tb.removerows(i)
        #for row in rowsToRemove:
            #tb.removerows(row)
        tb.close()

    # for Set if intent not in state table
    # then add a new row to the state table and change index (STATE_ID) in main table
    elif mode.lower() == 'set':
        # Keep track of the new value to set the state_id to
        newState = -1
        # Adding to intents col
        statetb = vis+'/STATE'
        tb.open(statetb, nomodify=False)
        intents = tb.getcol('OBS_MODE')
        # Check if the intent already exists
        if intent in intents:
            print("Intent already exists")
            newState = np.where(intents == intent)
        # If it doesn't add a row with the new intent
        else:
            tb.addrows(1)
            intents = list(tb.getcol('OBS_MODE'))
            intents[-1] = intent
            intents = np.asarray(intents)
            tb.putcol('OBS_MODE', intents)
            newState = len(intents) - 1
            tb.close()
            
            # For all selected rows replace with new state_id
            tb.open(vis, nomodify=False)
            stateCol = tb.getcol('STATE_ID')
            for row in selectedRows:
                stateCol[row] = newState
            tb.putcol('STATE_ID', stateCol)
            tb.close()
            
        tb.close()
    
    # For Append mode
    elif mode.lower() == 'append':
        statetb = vis+'/STATE'
        # Find our selected intents
        for i in selectedIntents:
            newState = -1
            tb.open(statetb, nomodify=False)
            intents = tb.getcol('OBS_MODE')
            # Add a row with old intent + new
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
        tb.open(vis, nomodify=False)
        stateCol = tb.getcol('STATE_ID')
        for row in selectedRows:
            if stateCol[row] in selectedIntents:
                stateCol[row] = selectedIntents[stateCol[row]]
        tb.putcol('STATE_ID', stateCol)
        tb.close()

    return changeList
