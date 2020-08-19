'''
MultiLoop

Member variable -- data type -- description:
self._parentLabel -- str -- parent label for all the multiloop subcomponents
self._subunitLabels -- list -- list of all the subunit labels for the multiloop
self._numSubunits -- int -- number of subunits composing the multiloop
self._sequences -- dictionary -- dictionary of the multiloop component sequences. key values are the subunit labels
self._span -- dictionary -- dictionary of the multiloop component spans. key values are the subunit labels
self._closingPairs -- ((str, str), (str, str)) -- tuple containing the 5' and 3' closing base pairs as tuples
self._closingPairsSpan - ((int, int), (int, int)) -- tuple containing the 5' and 3' closing base pair spans as tuples
'''


class Multiloop:
    # __init__() method for MultiLoop class
    def __init__(self, parentLabel, subunitLabels, sequences, spans, closingPairs, closingPairsSpan):
        self._parentLabel = parentLabel
        self._subunitLabels = subunitLabels
        self._numSubunits = len(sequences)
        self._sequences = sequences
        self._spans = spans
        self._closingPairs = closingPairs
        self._closingPairsSpan = closingPairsSpan

    ###
    # Internal Methods
    ###

    # define string representation for MultiLoop object
    def __str__(self):
        return f'MultiLoop: {self._parentLabel}'

    # internal method to set the 5' and 3' neighbors for a multiloop
    def _addNeighbors(self, neighbor5p, neighbor3p):
        self._neighbor5p = neighbor5p
        self._neighbor3p = neighbor3p

    ###
    # User Accesible Methods
    ###

    '''
    Function: MultiLoop.label()
    Description: Function to return the parent Label for the MultiLoop object
    Parameters:
            (newLabel=None) -- str -- new label to identify the Multiloop object
    Return Value:
            str - the label for the current MultlLoop object
    '''

    def label(self, newLabel=None):
        if newLabel:
            self._parentLabel = newLabel
        else:
            return self._parentLabel

    '''
    Function: MultiLoop.subunitLabels()
    Description: Function to return a list of the MultiLoop object subcomponents
    Parameters: None
    Return Value:
            list - list of the subunit labels for the MultiLoop object
    '''

    def subunitLabels(self):
        return self._subunitLabels

    '''
    Function: MultiLoop.numSubunits()
    Description: Function to return the number of subcomponents composing the MultiLoop
    Parameters: None
    Return Value:
            int - the number of subunits that compose the MultiLoop Structure
    '''

    def numSubunits(self):
        return self._numSubunits

    '''
    Function: MultiLoop.sequence()
    Description: Function to return dictionary of subunitLabel : Sequence pairs for the MultiLoop object
    Parameters:
            (subunit=None) -- str -- label for the specific subunit being accessed.
    Return Value:
            dict - dictionary of multiloop sequences mapped to their subunit label
            * if a specific subunit label is provided, only that sequence will bw returned
    '''

    def sequence(self, subunit=None):
        if(subunit):
            try:
                sequence = self._sequences[subunit]
                return sequence
            except KeyError:
                return None
        else:
            return self._sequences

    '''
    Function: MultiLoop.span()
    Description: Function to return dictionary of subunitLabel : SequenceSpans pairs for the MultiLoop object
    Parameters:
            (subunit=None) -- str -- label for specific subunit being accessed
    Return Value:
            dict - dictionary of tuples defining start and stop indicies for multiloop subunits mapped to their subunit label
            * if a specific subunit label is provided, only that tuple will be returned
    '''

    def span(self, subunit=None):
        if(subunit):
            try:
                span = self._spans[subunit]
                return span
            except KeyError:
                return None
        else:
            return self._spans

    '''
    Function: MultiLoop.closingPairs()
    Description: Function to return dictionary of subunitLabel : closingPairs for the MultiLoop object
    Parameters:
            (subunit=None) -- str -- label for specific subunit being accessed
    Return Value:
            dict - dictionary of tuples defining closing base pairs for multiloop subunits mapped to their subunit label
            * if a specific subunit label is provided, only that tuple will be returned
    '''

    def closingPairs(self, subunit=None):
        if(subunit):
            try:
                closingPair = self._closingPairs[subunit]
                return closingPair
            except KeyError:
                return None
        else:
            return self._closingPairs

    '''
    Function: MultiLoop.closingPairsSpan()
    Description: Function to return dictionary of subunitLabel : closingPairsSpan for the MultiLoop object
    Parameters:
            (subunit=None) -- str -- label for specific subunit being accessed
    Return Value:
            dict - dictionary of tuples defining closing base pair index locations for multiloop subunits mapped to their subunit label
            * if a specific subunit label is provided, only that tuple will be returned
    '''

    def closingPairsSpan(self, subunit=None):
        if(subunit):
            try:
                closingPairsSpan = self._closingPairsSpan[subunit]
                return closingPairsSpan
            except KeyError:
                return None
        else:
            return self._closingPairsSpan
