import numpy as np

## Free Energy Parameter Imports ##
# initiation parameters for internal loops, bulges, and hairpins
from bpRNAStructure.parameters.LoopInitiationEnergy import BulgeInit
# Watson-Crick stacking interaction parameters
from bpRNAStructure.parameters.StackingEnergies import StackingEnergies

## Free Energy Parameter Constants ##
R = 0.001987204258  # source: https://en.wikipedia.org/wiki/Gas_constant
T = 310.15

# Stems(source: https://rna.urmc.rochester.edu/NNDB/turner04/wc-parameters.html)
INTERMOLECULAR_INIT = 4.09  # intermolecular initiation value
STEM_SYMMETRY_PENALTY = 0.43
STEM_AU_END_PENALTY = 0.45

# Bulges
SPECIAL_C_BULGE = -0.9
BULGE_AU_END_PENALTY = 0.45

# other Constants
CANONICAL_BASE_PAIRS = [('A', 'U'), ('U', 'A'), ('G', 'C'),
                        ('C', 'G'), ('G', 'U'), ('U', 'G')]
'''
BULGE OBJECT
The Bulge object is used to represent the bulge RNA secondary structure.

Member Variable -- Data Type -- Description:
self._label -- string -- The label for the bulge as defined by the structure
    type file.
self._sequence -- string -- The RNA sequence for the bulge.
self._sequenceLen -- Int -- The length of the bulge as measured in nucleotides.
self._span -- (int, int) -- Tuple containing the integer start and stop
    indices for the RNA sequene that defines the bulge.
self._closingPair5p -- (string, string) -- Tuple containing 2 single
    character strings. The first string the the 5' base in 5' closing pair for
    the bule. The second character is the 3' base in the 5' closing pair.
self._closingPair5pSpan -- (int, int) -- Tuple containing 2 integers. The first
    integer is the index of the 5' base in 5' closing pair for the bulge.
    The second integer is the index of the 3' base in the 5' closing pair
self._closingPair3p -- (string, string) -- Tuple containing 2 single character
    strings. The first string the the 5' base in 3' closing pair for the bule.
    The second character is the 3' base in the 3' closing pair.
self._closingPair3pSpan -- (int, int) -- Tuple containing 2 integers. The first
    integer is the index of the 5' base in 3' closing pair for the bule.
    The second integer is the index of the 3' base in the 3' closing pair
self._pk -- int -- the pseudoknot the bulge is a part of, if any(default
    value is None)
self._neighbot5p -- str -- label for the 5'neighbor of the bulge
self._neighbot3p -- str -- label for the 3'neighbor of the bulge



                    C
            5'   AGC UAG   3'
                 ||| |||
            3'   UCG-AUC   5'
                     ^ 3' closing pair = ('U', 'A')
                   ^
                    5' closing pair = ('C', 'G')


'''


class Bulge:
    # __init__ method for bulge object
    def __init__(self, label=None, sequence='', sequenceSpan=(-1, -1),
                 closingPair5p=('', ''), closingPair5pSpan=(-1, -1),
                 closingPair3p=('', ''), closingPair3pSpan=(-1, -1),
                 pk=None, neighbor5p=None, neighbor3p=None):
        self._label = label
        self._sequence = sequence
        self._sequenceLen = len(sequence)
        self._span = sequenceSpan
        self._closingPair5p = closingPair5p
        self._closingPair5pSpan = closingPair5pSpan
        self._closingPair3p = closingPair3p
        self._closingPair3pSpan = closingPair3pSpan
        self._pk = pk
        self._neighbor5p = neighbor5p
        self._neighbor3p = neighbor3p

    ###
    # Internal Methods
    ###

    # defines string representation for object
    def __str__(self):
        return f'Bulge: {self._label}'

    # define len function operation for Bulge objects
    def __len__(self):
        return self._sequenceLen

    # internal method to set the 5' and 3' neighbors for a bulge
    def _addNeighbors(self, neighbor5p, neighbor3p):
        self._neighbor5p = neighbor5p
        self._neighbor3p = neighbor3p

    ###
    # User Accesible Functions
    ###

    '''
    Function: Bulge.label()
    Description: Function returns the label for the bulge object. Also allows
                user to define new label
    Parameters:
    Return Value:
    '''

    def label(self, newLabel=None):
        if newLabel:
            self._label = newLabel
        else:
            return self._label

    '''
    Function: Bulge.sequence()
    Description: Function returns the sequence that defines the bulge
                structure. Also allows user to define new sequence
    Parameters:
            (newSequence=None) -- str -- str nucleotides sequence to define
             new bulge sequence
    Return Value:
            str - the current nucleotide sequence that defines the bulge
    '''

    def sequence(self, newSequence=None):
        if newSequence:
            self._sequence = newSequence
            self._sequenceLen = len(newSequence)
        else:
            return self._sequence

    '''
    Function: Bulge.span()
    Description: Function returns the start and stop indices of the bulge as a
                tuple. Ex: (start, stop)
    Parameters: None
    Return Value:
            (int, int) - tuple containing the start and stop index positions
            of the bulge loop
    '''

    def span(self):
        return self._span

    '''
    Function: Bulge.sequenceLen()
    Description: Function returns the length of the bulge
    Parameters: None
    Return Value:
            int - the integer value length of the bulge
    '''

    def length(self):
        return self._sequenceLen

    '''
    Function: Bulge.closingPair5p()
    Description: Function returns a tuple containg the 5' closing pair for the
                 bulge. Also allows user to define new closing pair
    Parameters:
            (newClose=None) -- (str, str) -- new 5' closing pair for the
            bulge object
    Return Value:
            (str, str) - a tuple containing the 5' closing base pairs for
            the bulge object
    '''

    def closingPair5p(self, newClose=None):
        if newClose:
            self._closingPair5p = newClose
        else:
            return self._closingPair5p

    '''
    Function: Bulge.closingPair5pSpan()
    Description: Function returns a tuple containg the indices of the 5'
                closing pair for the bulge
    Parameters: None
    Return Value:
            (int, int) - tuple containing the index locations of the 5'
            closing base pair
    '''

    def closingPair5pSpan(self):
        return self._closingPair5pSpan

    '''
    Function: Bulge.closingPair3p(newClose=None)
    Description: Function returns a tuple containg the 3' closing pair for
                the bulge. Also allows user to define new closing pair
    Parameters:
            (newClose=None) -- (str, str) -- new 3' closing pair for the
            bulge object
    Return Value:
            (str, str) - a tuple containing the 3' closing base pairs for
            the bulge object
    '''

    def closingPair3p(self, newClose=None):
        if newClose:
            self._closingPair3p = newClose
        return self._closingPair3p

    '''
    Function: Bulge.closingPair5pSpan()
    Description: Function returns a tuple containg the indices of the 5'
                closing pair for the bulge
    Parameters: None
    Return Value:
            (int, int) - tuple containing the index locations of the 5'
            closing base pair
    '''

    def closingPair3pSpan(self):
        return self._closingPair3pSpan

    '''
    Function: Bulge.neighbors()
    Description: function to get the StructureComponents directly adjacent
                to the bulge
    Parameters: None
    Return Value:
            (str, str) - tuple containg the labels for the 5' and 3'
            neighbors of the bulge
    '''

    def neighbors(self):
        return (self._neighbor5p, self._neighbor3p)

    '''
    Function: Bulge.canonical()
    Description: function to check for valid conditions for calculating
                bulge energy
    Parameters: None
    Return Value:
            bool - returns True if all valid energy parameters are present
            for energy calculation
    '''

    def canonical(self):
        if self._sequenceLen == 1:
            if (self._closingPair5p not in StackingEnergies) or (self._closingPair3p not in StackingEnergies[self._closingPair5p]):
                return False
        return True

    '''
    Function: Bulge.energy()
    Description: function calculates the folding free energy change for the
                bulge
    Parameters:
            (strict=True) -- bool -- when true only energy values for bulges
            with all valid energy parameters will be calaculated
    Return Value:
            float - the calculated energy of the Bulge
    '''

    def energy(self, strict=True):
        if self._sequenceLen == 1:  # bulges of length 1
            # get base pair stack
            try:
                basePairStack = StackingEnergies[self._closingPair5p][self._closingPair3p]
            except KeyError:
                if strict:
                    return None  # strict mode - only calculate energy for bulges with valid params
                else:
                    basePairStack = 0

            # check for special C bulge case
            specialC = 0
            cCount = 0  # number of adjacent C's
            if self._sequence == 'C' and (self._closingPair5p[0] == 'C' or self._closingPair3p[0] == 'C'):
                specialC = SPECIAL_C_BULGE
                cCount = 1  # number of possible states due to adjacent C's
                if(self._closingPair5p[0] == 'C'):
                    cCount += 1
                if (self._closingPair3p[0] == 'C'):
                    cCount += 1

                return BulgeInit[1] + basePairStack + specialC - (R * T * np.log(cCount))

            # if not special C bulge, return bulge init + basePairStack
            else:
                return BulgeInit[1] + basePairStack

        else:  # bulge of length > 1
            if self._sequenceLen in BulgeInit:  # try to get value from dictionary
                return BulgeInit[self._sequenceLen]
            else:  # otherwise calculate
                return BulgeInit[6] + (1.75 * R * T * np.log(float(self._sequenceLen/6.0)))
