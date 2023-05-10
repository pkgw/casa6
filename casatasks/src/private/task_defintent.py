from __future__ import absolute_import
import os
import numpy as np

import casatools
    
def selectionToQuery(queryStrings):
    # querys are constructed from a list of the query parameters
    fieldQuery = queryStrings[0]
    scanQuery = queryStrings[1]
    obsQuery = queryStrings[2]
    
    queryTypes = ['FIELD_ID', 'SCAN_NUMBER', 'OBSERVATION_ID']
    
    querys = []
    finalQuery = ''
    
    for i in range(len(queryTypes)):
        queryList = []
        # split the query string for this selection
        elements = queryStrings[i].split(',')
        
        for element in elements:
            if '~' in element:
                start, end = element.split('~')
                start = int(start.strip())
                end = int(end.strip())
                
                queryList += [x for x in range(start, end+1)]
            elif element != '':
                queryList.append(int(element.strip()))
        
        if queryList != []:
            querys.append(f"{queryTypes[i]} in {queryList}")
    if querys != []:
        finalQuery = " && ".join(querys)
        
    return finalQuery
            

def defintent(vis='', intent='', mode='',
              scan='', field='', obsid=''):
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
    """

    tb = casatools.table()
    
    # If no intent has been provided exit the task and print
    if vis == '':
        print('You must specify a MS')
        return
        
    if field == '':
        print('You must specify a field ID or name')
        
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
    tmpField = list(np.where(names==fieldSplit)[0])
    tmpField = [str(x) for x in tmpField]
    field = ",".join(tmpField)
            
    # Get exitsing intents
    tb.open(vis+'/STATE')
    intentcol = tb.getcol('OBS_MODE')
    tb.close()
    
    # Allowed Intents? / Intents reformatting?
    
    # If selected field is found
    foundField = False
    if (type(scan) != list):
        scan = str(scan)
    
    selectedRows = set()
    selectedIntents = dict()
    
    tb.open(vis)
    fieldIds = tb.getcol('FIELD_ID')
    scanNum = tb.getcol('SCAN_NUMBER')
    stateIds = tb.getcol('STATE_ID')
    obsIds = tb.getcol('OBSERVATION_ID')
    
    # query tool selection
    taskQuery = selectionToQuery([field, scan, obsid])
    
    selectedData = tb.query(taskQuery)
    selectedRows = set(selectedData.rownumbers())
    
    selectedStateIds = selectedData.getcol('STATE_ID')
    for i in range(len(selectedRows)):
        selectedIntents[selectedStateIds[i]] = selectedStateIds[i]
        
    tb.close()
                
    print("Number of matching rows found: ", len(selectedRows))
    
    # for Set if intent not in state table
    # then add a new row to the state table and change index (STATE_ID) in main table
    if mode.lower() == 'set':
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
            intents = tb.getcol('OBS_MODE')
            intents[-1] = intent
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
    if mode.lower() == 'append':
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
                intents = tb.getcol('OBS_MODE')
                intents[-1] = newIntent
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

    return
