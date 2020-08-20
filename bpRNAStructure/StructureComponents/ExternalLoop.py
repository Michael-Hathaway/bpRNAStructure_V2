'''
External Loops

Member variable -- data type -- description:
self._label -- string -- label for the external loop secondary structure
self._sequence -- string -- base sequence that defines the external loop
self._sequenceLen -- int -- length of the external loop sequence
self._span --tuple(int, int) -- tuple containing the integer start and stop
    locations for the external loop sequence
self._closingPair5p -- tuple(string, string) -- tuple that contains the 5'
    closing base pair for the external loop
self._closingPair5pSpan -- tuple(int, int) -- tuple containg the integer
    index locations for the 5' closing pair
self._closingPair3p -- tuple(string, string) -- tuple that contains the 3'
    closing base pair for the external loop
self._closingPair3pSpan -- tuple(int, int) -- tuple containg the integer
    index locations for the 3' closing pair
'''


class ExternalLoop:
    # __init__() method for the external loop object
    def __init__(self, label='', seq='', seqSpan=(-1, -1),
                 closingPair5p=('', ''), closingPair5pSpan=(-1, -1),
                 closingPair3p=('', ''), closingPair3pSpan=(-1, -1),
                 neighbor5p=None, neighbor3p=None):
        self._label = label
        self._sequence = seq
        self._sequenceLen = len(seq)
        self._span = seqSpan
        self._closingPair5p = closingPair5p
        self._closingPair5pSpan = closingPair5pSpan
        self._closingPair3p = closingPair3p
        self._closingPair3pSpan = closingPair3pSpan
        self._neighbor5p = neighbor5p
        self._neighbor3p = neighbor3p

    ###
    # Internal Methods
    ###

    # Defines the string representation of the external loop
    def __str__(self):
        return f'External Loop: {self._label}'

    # define len function operation for ExternalLoop objects
    def __len__(self):
        return self._sequenceLen

    # internal method to set the 5' and 3' neighbors for a external loop
    def _addNeighbors(self, neighbor5p, neighbor3p):
        self._neighbor5p = neighbor5p
        self._neighbor3p = neighbor3p

    ###
    # User Accesible Methods
    ###

    '''
    Function: ExternalLoop.label()
    Description: Function returns the label for the external loop
    Parameters:
            (newLabel=None) -- str -- new label to define the ExternalLoop
            Object
    Return Value:
            str - the current label for the ExternalLoop object
    '''

    def label(self, newLabel=None):
        if newLabel:
            self._label = newLabel
        else:
            return self._label

    '''
    Function: ExternalLoop.sequence()
    Description: Function returns the sequence that defines the external loop
    Parameters:
            (newSequence=None) -- str -- new nucleotide sequence to define the
            ExternalLoop
    Return Value:
            str - current sequence that defines the ExternalLoop
    '''

    def sequence(self):
        return self._sequence

    '''
    Function: ExternalLoop.sequenceLen()
    Description: function returns the length of the external loop sequence
    Parameters: None
    Return Value:
            int - the integer value length of the ExternalLoop sequence
    '''

    def sequenceLen(self):
        return self._sequenceLen

    '''
    Function: ExternalLoop.span()
    Description: Function returns a tuple containing the start and stop index
        locations for the external loop sequence
    Parameters: None
    Return Value:
            (int, int) - tuple containing the start and stop index locations
            of the ExternalLoop
    '''

    def span(self):
        return self._span

    '''
    Function: ExternalLoop.neighbors()
    Description: function to get the StructureComponents directly adjacent to
        the external loop
    Parameters: None
    Return Value:
            (str, str) - tuple containing the labels for the 5' and 3'
            neighbors of the ExternalLoop
    '''

    def neighbors(self):
        return (self._neighbor5p, self._neighbor3p)
