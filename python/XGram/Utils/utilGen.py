'''Generic utility functions.'''

def join(pDict):
    '''Merges the string lists of a dictionary.

    Input: 
       pDict = a dictionary
    Return:
       Mutated dictionary'''
    for key in pDict:
	pDict[key] = "".join(pDict[key])

def append(pDict, pKey, pList):
    '''Appends list to a dictionary item, if necessary creating the list.

    Input: 
       pDict = a dictionary
       pKey = key value for item
       pList = list to be appended
    Return:
       Mutated dictionary'''
    if pKey not in pDict: pDict[pKey] = []
    pDict[pKey] += pList

