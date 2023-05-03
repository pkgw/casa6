from __future__ import absolute_import
import os
import numpy as np

import casatools

def parseSelection(selection):
    selected = []
    
    if selection == '':
        return selected
        
    for item in selection.split(','):
        if '~' in item:
            start, end = item.split('~')
            selected.extend(range(int(start), int(end)+1))
        else:
            selected.append(int(item))
    return selected
    
def selectionToQuery(queryParam, queryString):
    elements = queryString.split(',')
    
    queryElements = []
    for element in elements:
        if '~' in element:
            start, end = element.split('~')
            start = int(start.strip())
            end = int(end.strip())
            
            rangeString = f"{start} <= {queryParam} <= {end}"
            queryElements.append(rangeString)
        else:
            queryElements.append(f"{queryParam} == {element.strip()}")
            
    queryFinal = " || ".join(queryElements)
    
    return queryFinal
            

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
    """
    # query selection
    toJoin = []
    if (field != ''):
        toJoin.append(selectionToQuery('FIELD_ID', field))
    if (scan != ''):
        toJoin.append(selectionToQuery('SCAN_NUMBER', scan))
    if (obsid != ''):
        toJoin.append(selectionToQuery('OBSERVATION_ID', obsid))
        
    print(toJoin)
    queryString = " && ".join(toJoin)
    selectedData = tb.query(queryString)
    
    selectedRows = set(selectedData.rownumbers())
    selectedStateIds = selectedData.getcol('STATE_ID')
    for row in selectedRows:
        selectedIntents[selectedStateIds[row]] = selectedStateIds[row]
    """
    tb.close()
    
    # split selection parameters into array
    selectedFieldList = parseSelection(field)
    selectedScanList = parseSelection(scan)
    selectedObsIdList = parseSelection(obsid)
    
    # Select rows based on field and scan and add selected intents
    # if row field/scan is in the sting array, select that row
    for row in range(len(fieldIds)):
        # also select if field == ''
        if field == '' or fieldIds[row] in selectedFieldList or fieldIds[row] in np.where(fieldnames == field):
            foundField = True
            if scanNum[row] in selectedScanList or scan == '':
                if obsIds[row] in selectedObsIdList or obsid == '':
                    #selectedRows[row] = stateIds[row]
                    selectedRows.add(row)
                    selectedIntents[stateIds[row]] = stateIds[row]
                    #selectedIntents.add(stateIds[row])
                
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
