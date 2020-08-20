# Free Energy Parameter Imports ##
# Watson-Crick stacking interaction parameters
from bpRNAStructure.parameters.StackingEnergies import StackingEnergies
# stacking terminal mismatches for Hairpin calculations
from bpRNAStructure.parameters.StackTerminalMismatches import StackTerminalMismatches

# Free Energy Parameter Constants ##
R = 0.001987204258  # source: https://en.wikipedia.org/wiki/Gas_constant
T = 310.15

# Stems(source: https://rna.urmc.rochester.edu/NNDB/turner04/wc-parameters.html)
INTERMOLECULAR_INIT = 4.09  # intermolecular initiation value
STEM_SYMMETRY_PENALTY = 0.43
STEM_AU_END_PENALTY = 0.45

# other Constants
CANONICAL_BASE_PAIRS = [('A', 'U'), ('U', 'A'), ('G', 'C'),
                        ('C', 'G'), ('G', 'U'), ('U', 'G')]
'''
## STEM OBJECT ##
the Stem object is used to represent RNA secondary structure stems.

Member variable -- data type -- description:
self._label -- String -- the label for the stem as defined in the structure
    type file.
self._sequence5p -- String -- the 5' portion of the stem sequence.
self._sequence3p -- String -- the 3' portion of the stem sequence.
self._sequenceLen -- Int -- the length of the stem in number of base pairs.
self._sequence5p_index -- (int, int) -- tuple containing the integer value
    start and stop indices for the 5' portion of the stem sequence.
self._sequence3p_index -- (int, int) -- tuple containing the integer value
    start and stop indices for the 3' portion of the stem sequence.
self._neighbor5p -- str -- label for 5' neighbor in Structure object
self._neighbor3p -- str -- label for 5' neighbor in Structure object


            5' Sequence

            5' ACGUG 3'
               |||||
            3' UGCAC 5'

            3' Sequence

'''


class Stem:
    # __init__ method for stem object
    def __init__(self, label="", sequence5p="", sequence3p="",
                 sequence5pSpan=(-1, -1), sequence3pSpan=(-1, -1),
                 neighbor5p=('', ''), neighbor3p=('', ''),
                 adjacentBulges=(False, False)):
        self._label = label  # sequence label
        self._sequence5p = sequence5p  # 5' portion of stem
        self._sequence3p = sequence3p  # 3' portion of stem
        self._sequence = list(
            zip(list(self._sequence5p), list(self._sequence3p[::-1])))
        self._sequenceLen = (
            len(sequence5p) + len(sequence3p)) // 2  # sequence length
        # tuple containing start and stop indices of 5' prime portion of stem
        self._sequence5pSpan = sequence5pSpan
        # tuple containing start and stop indices of 3' prime portion of stem
        self._sequence3pSpan = sequence3pSpan
        self._neighbor5p = neighbor5p
        self._neighbor3p = neighbor3p
        self._adjacentBulges = adjacentBulges

    ###
    # Internal Methods
    ###

    # define string representation of object
    def __str__(self):
        return f'Stem: {self._label}'

    # define len function operation for Stem objects
    def __len__(self):
        return self._sequenceLen

    # Internal method to set the object _sequence member variable as a list of
    # tuples based on the 5' and 3' sequence variables
    def _setSequence(self):
        if len(self._sequence5p) == len(self._sequence3p):
            self._sequence = list(
                zip(list(self._sequence5p), list(self._sequence3p[::-1])))

    # internal method to update the sequenceLen member variable when the
    # sequence is changed by the user
    def _setSequenceLen(self):
        self._sequenceLen = (len(self._sequence5p) +
                             len(self._sequence3p)) // 2

    # internal method used during Structure object parsing to track if
    # stem is next to length=1 bulges
    def _addAdjacentBulgeBoolean(self, bulge5p, bulge3p):
        self._adjacentBulges = (bulge5p, bulge3p)

    # Internal method that returns tuple containg booleans for whether or
    # not the stem is adjacent to length=1 bulges
    def _adjacentBulgeBoolean(self):
        return self._adjacentBulges

    # internal method to set the 5' and 3' neighbors for a stem
    def _addNeighbors(self, neighbor5p, neighbor3p):
        self._neighbor5p = neighbor5p
        self._neighbor3p = neighbor3p

    ###
    # User Accesible Methods
    ###

    '''
    Function: Stem.label()
    Description: function returns the label for the stem object. Also allows
                for user to change label of stems
    Parameters:
            (newLabel=None) -- str -- new label that to update the stem._label
            member variable
    Return Value:
            str - label for the Stem object
    '''

    def label(self, newLabel=None):
        if newLabel:
            self._label = newLabel
        else:
            return self._label

    '''
    Function: Stem.sequence5p()
    Description: function returns the 5' portion of the stem sequence
    Parameters:
            (newSequence) -- str -- new RNA sequence to define the 5' sequence
            of the stem.
    Return Value:
            str - The current 5' sequence for the Stem object
    '''

    def sequence5p(self, newSequence=None):
        if (newSequence):  # if new sequence is provided
            # check that sequence length matchees other 3' sequences
            if(len(newSequence) == self._sequenceLen):
                self._sequence5p = newSequence
                self._setSequence()  # reset the self._sequence variable
            else:
                print('Unable to set new 5\' sequence')
        else:
            return self._sequence5p

    '''
    Function: Stem.sequence5p()
    Description: function returns the 3' portion of the stem sequence
    Parameters:
            (newSequence) -- str -- new RNA sequence to define the 3'
            sequence of the stem.
    Return Value:
            str - The current 3' sequence for the Stem object
    '''

    def sequence3p(self, newSequence=None):
        if (newSequence):  # if new sequence is provided
            # check that sequence length matchees other 5' sequences
            if(len(newSequence) == self._sequenceLen):
                self._sequence3p = newSequence
                self._setSequence()  # reset the self._sequence variable
            else:
                print('Unable to set new 3\' sequence')
        else:
            return self._sequence3p

    '''
    Function: Stem.sequence()
    Description: function returns the stem sequence as a list of tuples
                containg base pairs. Ex: [('C','G'), ... , ('A', 'U')]
    Parameters:
            (sequence5p=None) -- str -- new RNA sequence to define the 5'
            sequence of the stem.
            (sequence3p=None) -- str -- new RNA sequence to define the 3'
            sequence of the stem.
    Return Value:
            list - list of tuples representing the base pair sequence of
            the stem.
    '''

    def sequence(self, sequence5p=None, sequence3p=None):
        if(sequence5p and sequence3p):
            if(len(sequence5p) == len(sequence3p)):
                self._sequence5p = sequence5p
                self._sequence3p = sequence3p
                self._setSequence()
                self._setSequenceLen()
            else:
                print(
                    'Invalid - sequences are different lengths.')
        else:
            return self._sequence

    '''
    Function: Stem.sequenceLen()
    Description: Function returns the length of the Stem object
    Parameters: None
    Return Value:
            int - the integer value length of the stem
    '''

    def sequenceLen(self):
        return self._sequenceLen

    '''
    Function: Stem.span()
    Description: function returns a tuple containing two tuples that contain
                start and stop indices for the 5' and 3' sequence of the stem
    Parameters: None
    Return Value:
            ((int, int), (int, int)) - a tuple containing the tuple(int, int)
            start and stop positions for the 5' and 3' stem sequences
    '''

    def span(self):
        return (self._sequence5pSpan, self._sequence3pSpan)

    '''
    Function: Stem.sequence5pSpan()
    Description: function returns the start and stop indices of the 5' portion
                of the stem in a tuple. Ex: (start, stop)
    Parameters: None
    Return Value:
            (int, int) - a tuple containing the integer value start and stop
            positions of the 5' portion of the stem
    '''

    def sequence5pSpan(self):
        return self._sequence5pSpan

    '''
    Function: Stem.sequence3pSpan()
    Description: function returns the start and stop indices of the 3' portion
                of the stem in a tuple. Ex: (start, stop)
    Parameters: None
    Return Value:
            (int, int) - a tuple containing the integer value start and stop
            positions of the 3' portion of the stem
    '''

    def sequence3pSpan(self):
        return self._sequence3pSpan

    '''
    Function: Stem.neighbors()
    Description: Function returns a tuple containing the labels for the 5'
                and 3' neighbors of the stem
    Parameters: None
    Return Value:
            ((str, str), (str, str)) - tuple containging tuples containing
            the string value labels for the neighbors structures of the 5'
            and 3' portions of the stem sequennce
    '''

    def neighbors(self):
        return (self._neighbor5p, self._neighbor3p)

    '''
    Function: Stem.cannonical()
    Description: Function to check if all base pairs in a stem are canonical
                base pairings
    Parameters: None
    Return Value:
            bool - true or false as to whether or not the stem contains all
            cannonical base pairings
    '''

    def canonical(self):
        return (self._sequenceLen > 1 and all(pair in CANONICAL_BASE_PAIRS for pair in self._sequence))

    '''
    Function: Stem.energy()
    Description: function calculates the folding free energy change for the
                stem
    Parameters:
            (strict=True) -- bool -- when true, energy values will only be
            calculated for cannonical stems/stems with all present energy
            parameters
            (init=False) -- bool -- when true, the 4.09 Kcal/mol initiation
            value is inlcuded in energy calculations.
    Return Value:
            float - the calculated energy value for the given stem
    '''

    def energy(self, strict=True, init=False):
        if(self._sequenceLen == 1):
            return None

        seq = self.sequence()  # get stem as list of tuple base pairs

        # check for symmetry
        symmetry = 0
        if self._sequence5p == self._sequence3p:
            symmetry = STEM_SYMMETRY_PENALTY

        # check for AU end penalty
        endPenalty = 0
        if (seq[0] == ('A', 'U') or seq[0] == ('U', 'A') or seq[0] == ('G', 'U') or seq[0] == ('U', 'G')) and (self._adjacentBulgeBoolean()[0] == False):
            endPenalty += STEM_AU_END_PENALTY
        if (seq[-1] == ('A', 'U') or seq[-1] == ('U', 'A') or seq[-1] == ('G', 'U') or seq[-1] == ('U', 'G')) and (self._adjacentBulgeBoolean()[1] == False):
            endPenalty += STEM_AU_END_PENALTY

        # sum up watson crick stacking interactions
        stack = 0
        for i in range(0, self._sequenceLen-1):
            try:
                stack += StackingEnergies[seq[i]][seq[i+1]]
            except KeyError:
                if strict:  # default strict mode - only calculate energy for stems with all valid parameters
                    return None
                else:
                    continue

        if(init):
            return INTERMOLECULAR_INIT + symmetry + endPenalty + stack
        else:
            return symmetry + endPenalty + stack
