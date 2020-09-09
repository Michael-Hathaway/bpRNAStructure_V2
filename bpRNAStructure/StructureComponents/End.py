'''
ENDS

Member variable -- data type -- description:
self._label -- string -- label for the end objects
self._sequence -- string -- sequence that defines the end objects
self._span -- tuple(int, int) -- tuple containing the integer start and
    stop locations for the end object
'''


class End:
    # __init__() method for end object
    def __init__(self, label='', sequence='', span=(-1, -1), neighbor=None):
        self._label = label
        self._sequence = sequence
        self._sequenceLen = len(sequence)
        self._span = span
        self._neighbor = None

    ###
    # Internal Methods
    ###

    # define string representation of end object
    def __str__(self):
        return f'End: {self._label}'

    # define len function operation for End objects
    def __len__(self):
        return self._sequenceLen

    # internal method to set the 5' and 3' neighbors for a end
    def _addNeighbors(self, neighbor5p, neighbor3p):
        self._neighbor5p = neighbor5p
        self._neighbor3p = neighbor3p

    ###
    # User Accesible Methods
    ###

    '''
    Function: End.label()
    Description: Function returns the label for the end object
    Parameters:
            (newLabel=None) -- str -- new label to define the End object
    Return Value:
            str - the current label for the End object
    '''

    def label(self, newLabel=None):
        if newLabel:
            self._label = newLabel
        else:
            return self._label

    '''
    Function: End.sequence()
    Description: Function returns the sequence that defines the end object
    Parameters:
            (newSequence=None) -- str -- new nucleotide sequence to define
            the End object
    Return Value:
            str - the current sequence that defines the End object
    '''

    def sequence(self, newSequence=None):
        if newSequence:
            self._sequence = newSequence
        else:
            return self._sequence

    '''
    Function: End.sequenceLen()
    Description: function returns the length of the End sequence
    Parameters: None
    Return Value:
            int - the integer value length of the End object
    '''

    def length(self):
        return self._sequenceLen

    '''
    Function: End.span()
    Description: Function returns a tuple that contains the integer start and
        stop index locations for the end object
    Parameters: None
    Return Value:
            (int, int) - tuple containing the integer value start and stop
            indices of the End object
    '''

    def span(self):
        return self._span

    '''
    Function: End.neighbors()
    Description: function to get the StructureComponents directly adjacent to
        the end
    Parameters: None
    Return Value:
            (str, str) - tuple containing the labels for the
            StructureComponents neighboring the End in the Structure object
    '''

    def neighbors(self):
        return (self._neighbor5p, self._neighbor3p)
