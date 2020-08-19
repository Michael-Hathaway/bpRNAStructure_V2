'''
NON-CANONICAL BASE PAIRINGS

Member variable -- data type -- description:
self._label -- string -- label for the NCBP objects
self._basePair -- tuple(string, string) -- tuple containing the base pairs that define the NCBP object
self._basePairSpan -- tuple(int, int) -- tuple containing the integer locations of the NCBP
self._parentUnit -- string -- label for the secondary structure that the NCBP is located in
'''


class NCBP:
    # __init__() method for the NCBP object
    def __init__(self, label, basePair, basePairSpan, loc):
        self._label = label
        self._basePair = basePair
        self._basePairSpan = basePairSpan
        self._parentUnit = loc

    ###
    # Internal Methods
    ###

    # Defines the string representation of the NCBP object
    def __str__(self):
        return f'NCBP: {self._label}'

    ###
    # User Accesible Methods
    ###

    '''
    Function: NCBP.label()
    Description: Functions returns the label for the NCBP object
    Parameters:
            (newLabel=None) -- str -- new label to identify the NCBP object
    Return Value:
            str - the current label for the NCBP object
    '''

    def label(self, newLabel=None):
        if newLabel:
            self._label = newLabel
        else:
            return self._label

    '''
    Function: NCBP.pair()
    Description: Function returns a tuple containing the two base pairs that define the NCBP object
    Parameters:
            (newPair=None) -- (str, str) -- new tuple to define the NCBP
    Return Value:
            (str, str) - the current base pair that defines the NCBP
    '''

    def pair(self, newPair=None):
        if newPair:
            self._basePair = newPair
        else:
            return self._basePair

    '''
    Function: NCBP.span()
    Description: Function returns a tuple containing the integer locations of the base pairs that define the NCBP
    Parameters: None
    Return Value:
            (int, int) - tuple containing the index locations of the NCBP
    '''

    def span(self):
        return self._basePairSpan

    '''
    Function: NCBP.parentUnit()
    Description: Function returns a string the identifies the secondary structure that the NCBP occurs in
    Parameters: None
    Return Value:
            str - label for the StructureComponent that contains the NCBP
    '''

    def parentUnit(self):
        return self._parentUnit
