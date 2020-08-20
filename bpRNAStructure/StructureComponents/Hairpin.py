## Free Energy Parameter Imports ##
from bpRNAStructure.parameters.LoopInitiationEnergy import InternalLoopInit, BulgeInit, HairpinInit
from bpRNAStructure.parameters.StackingEnergies import StackingEnergies
from bpRNAStructure.parameters.InnerLoop_1x1_Energies import InnerLoop_1x1_Energies
from bpRNAStructure.parameters.InnerLoop_1x2_Energies import InnerLoop_1x2_Energies
from bpRNAStructure.parameters.InnerLoop_2x2_Energies import InnerLoop_2x2_Energies
from bpRNAStructure.parameters.InnerLoopMismatches import InnerLoopMismatches_2x3, OtherInnerLoopMismtaches
from bpRNAStructure.parameters.StackTerminalMismatches import StackTerminalMismatches
from bpRNAStructure.parameters.SpecialHairpins import SpecialHairpins

## Free Energy Parameter Constants ##
R = 0.001987204258  # source: https://en.wikipedia.org/wiki/Gas_constant
T = 310.15

# Stems(source: https://rna.urmc.rochester.edu/NNDB/turner04/wc-parameters.html)
INTERMOLECULAR_INIT = 4.09  # intermolecular initiation value
STEM_SYMMETRY_PENALTY = 0.43
STEM_AU_END_PENALTY = 0.45

# Inner loops
INNER_LOOP_ASYMMETRY_PENALTY = 0.6

# Bulges
SPECIAL_C_BULGE = -0.9
BULGE_AU_END_PENALTY = 0.45

# Hairpins
HAIRPIN_UU_GA_FIRST_MISMATCH_BONUS = -0.9
HAIRPIN_GG_FIRST_MISMATCH_BONUS = -0.8
HAIRPIN_SPECIAL_GU_CLOSURE = -2.2
HAIRPIN_C3 = 1.5
HAIRPIN_C_LOOP_A = 0.3
HAIRPIN_C_LOOP_B = 1.6

# other Constants
CANONICAL_BASE_PAIRS = [('A', 'U'), ('U', 'A'), ('G', 'C'),
                        ('C', 'G'), ('G', 'U'), ('U', 'G')]
'''
## HAIRPIN OBJECT ##
the Hairpin object is used to represent RNA secondary structure hairpins.

Member variable -- data type -- description:
self._label -- string -- the label for the hairpin as defined by the structure
    type file.
self._sequence -- string -- the RNA sequence for the hairpin.
self._sequenceLen -- Int -- the length of the hairpin as measured in number
    of nucleotides.
self._span -- (int, int) -- tuple containing the integer start and stop
    indices for the hairpin.
self._closingPair -- (string, string) -- tuple containing two single character
    strings. The first character corresponds to the 5' base in the closing
    pair. The second character is the 3' base in the closing pair.
self._closing_span -- (int, int) -- tuple containing two integers. The first
    integer is the index location of the 5' base in the closing pair. The second integer is the index location of the 3'base in the closing pair.
self._pk -- Int -- The pseudoknot the hairpin is a part of, if any(default
    value is None)
self._neighbor -- str -- label for the neighboring stem to the hairpin


                  C
                A   G
              G       A
               C     G <- first mismatch = ('C', 'G')
                A - U <- _Closing Pair = ('A', 'U')
                C - G
                G - C
                5'  3'


'''


class Hairpin:
    # __init__ method for stem object
    def __init__(self, label="", sequence="", sequenceSpan=(-1, -1),
                 closingPair=('', ''), closingPairSpan=(-1, -1),
                 pk=None, neighbors=('', '')):
        self._label = label
        self._sequence = sequence
        self._sequenceLen = len(sequence)
        self._span = sequenceSpan
        self._closingPair = closingPair
        self._closingPairSpan = closingPairSpan
        self._pk = pk
        self._neighbors = neighbors

    ###
    # Internal Methods
    ###

    # define string representation of object

    def __str__(self):
        return f'Hairpin: {self._label}'

    # define len function operation for Hairpin objects
    def __len__(self):
        return self._sequenceLen

    # Internal function to add hairpin neighbors to the object
    def _addNeighbors(self, neighbor5p, neighbor3p):
        self._neighbors = (neighbor5p, neighbor3p)

    ###
    # User Accessible Methods
    ###

    '''
    Function: Hairpin.label()
    Description: Function returns the label for the hairpin object. also
        allows user to define new label
    Parameters:
            (newLabel=None) -- str -- new label to indentify the hairpin object
    Return Value:
            str - the current label for the given hairpin object
    '''

    def label(self, newLabel=None):
        if newLabel:
            self._label = newLabel
        else:
            return self._label

    '''
    Function: Hairpin.sequence()
    Description: Function returns the sequence that defines the hairpin
        structure. Also allows user to define new sequence
    Parameters:
            (newSequence=None) -- str -- new sequence to define the hairpin
            loop
    Return Value:
            str - the current sequence that defines the hairpin loop
    '''

    def sequence(self, newSequence=None):
        if newSequence:
            self._sequence = newSequence  # set new sequence
            self._sequenceLen = len(newSequence)  # update sequence length
        else:
            return self._sequence

    '''
    Function: Hairpin.SequenceLen()
    Description: function returns the length of the hairpin
    Parameters: None
    Return Value:
            int - the length of the hairpin sequence
    '''

    def sequenceLen(self):
        return self._sequenceLen

    '''
    Function: Hairpin.span()
    Description: function returns the start and stop indices of the hairpin
        as a tuple. Ex: (start, stop)
    Parameters: None
    Return Value:
            (int, int) - tuple containing the integer value start and stop
            positions of the hairpin loop
    '''

    def span(self):
        return self._span

    '''
    Function: Hairpin.closingPair()
    Description: function returns a tuple that contains the closing pair for
        the hairpin. Ex: (5' closing base, 3' closing base). Also allows user
        to define new closing pair
    Parameters:
            (newClose) -- (str, str) -- new closing pair for the Hairpin object
    Return Value:
            (str, str) - tuple containing the 5' and 3' closing bases for
            the Hairpin
    '''

    def closingPair(self, newClose=None):
        if newClose:
            try:
                if newClose[0] and newClose[1]:
                    self._closingPair = newClose
            except:
                print(
                    'Invalid form provided for closing pair.')
        else:
            return self._closingPair

    '''
    Function: Hairpin.closingPairSpan()
    Description: Function returns the index locations of the closing pair
        bases as a tuple. Ex: (5' closing index, 3' closing index)
    Parameters: None
    Return Value:
            (int, int) - tuple containing the integer value locations of the
            closing base pairs for the Hairpin
    '''

    def closingPairSpan(self):
        return self._closingPairSpan

    '''
    Function: Hairpin.hairpinPK()
    Description: function returns the pseadoknot label for the hairpin if it
        exists
    Parameters: None
    Return Value:
            (int) - the pseudoknot that the hairpin is a part of if it exist.
            Will return none if not part of pseudoknot
    '''

    def hairpinPK(self):
        return self._pk

    '''
    Function: Hairpin.neighbors()
    Description: funtion to get the labels for the StructureComponents
                to the hairpin
    Parameters: None
    Return Value:
            (str, str) - tuple containing the labels for the neighboring
            secondary structures of the 5' and 3' ends of the sequence.
            Note: for hairpins the neighbors will always be the same.
    '''

    def neighbors(self):
        return self._neighbors

    '''
    Function: Hairpin.cannonical()
    Description: Function to check if the correct parameters are available
        to calculate the energy of the hairpin
    Parameters: None
    Return Value:
            bool - function return True if all necessary energy values are
            present for the Hairpin
    '''

    def canonical(self):
        firstMismatch = (self._sequence[0], self._sequence[-1])
        if(self._closingPair not in StackTerminalMismatches) or (firstMismatch not in StackTerminalMismatches[self._closingPair]):
            return False
        elif self._sequenceLen < 3:
            return False
        else:
            return True

    '''
    Function: Hairpin.energy()
    Description: function to calculate folding free energy of hairpin
    Parameters:
            (strict=True) -- bool -- when True, the function will only
            calculate the energy of the molecule valid energy parameters are
            present.
    Return Value:
            float - the calculated energy for the hairpin
    '''

    def energy(self, strict=True):
        # check that hairpin is at least 3 nucleotides long
        if self._sequenceLen < 3:
            return None

        # Check if the hairpin is a special case hairpin with precalculated energy values
        elif (self._closingPair in SpecialHairpins) and (self._sequence in SpecialHairpins[self._closingPair]):
            return SpecialHairpins[self._closingPair][self._sequence]

        # Hairpins of length 3
        elif self._sequenceLen == 3:
            # get hairpin initiation term
            if self._sequenceLen in HairpinInit:  # try to get from dictionary
                init = HairpinInit[self._sequenceLen]
            else:  # otherwise calculate
                init = HairpinInit[9] + (1.75 * R * T *
                                         np.log(float(self._sequenceLen/9.0)))

            # check for all c loop penalty
            if self._sequence.count('C') == self._sequenceLen:
                return init + HAIRPIN_C3

            return init

        # hairpins of 4 nucleotides or greater
        else:
            # get hairpin initiation term
            if self._sequenceLen in HairpinInit:  # try to get from dictionary
                init = HairpinInit[self._sequenceLen]
            else:  # otherwise calculate
                init = HairpinInit[9] + (1.75 * R * T *
                                         np.log(float(self._sequenceLen/9.0)))

            # get terminal mismatch parameter
            firstMismatch = (self._sequence[0], self._sequence[-1])
            try:
                terminalMismatch = StackTerminalMismatches[self._closingPair][firstMismatch]
            except KeyError:
                if strict:
                    return None  # strict mode - only calculate energy for hairpins with valid params
                else:
                    terminalMismatch = 0

            # UU/GA first mismatch bonus
            uu_ga_bonus = 0
            if firstMismatch == ('U', 'U') or firstMismatch == ('G', 'A'):
                uu_ga_bonus = HAIRPIN_UU_GA_FIRST_MISMATCH_BONUS

            # GG first mismatch
            gg_bonus = 0
            if firstMismatch == ('G', 'G'):
                gg_bonus = HAIRPIN_GG_FIRST_MISMATCH_BONUS

            # Special GU closure
            gu_closure = 0
            if self._closingPair == ('G', 'U') and firstMismatch == ('G', 'G'):
                gu_closure = HAIRPIN_SPECIAL_GU_CLOSURE

            # All C loop penalty
            c_loop_penalty = 0
            if self._sequence.count('C') == self._sequenceLen:
                c_loop_penalty = (self._sequenceLen *
                                  HAIRPIN_C_LOOP_A) + HAIRPIN_C_LOOP_B

            return init + terminalMismatch + uu_ga_bonus + gg_bonus + gu_closure + c_loop_penalty
