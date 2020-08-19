import numpy as np

## Free Energy Parameter Imports ##
from bpRNAStructure.parameters.LoopInitiationEnergy import InternalLoopInit
from bpRNAStructure.parameters.StackingEnergies import StackingEnergies
from bpRNAStructure.parameters.InnerLoop_1x1_Energies import InnerLoop_1x1_Energies
from bpRNAStructure.parameters.InnerLoop_1x2_Energies import InnerLoop_1x2_Energies
from bpRNAStructure.parameters.InnerLoop_2x2_Energies import InnerLoop_2x2_Energies
from bpRNAStructure.parameters.InnerLoopMismatches import InnerLoopMismatches_2x3, OtherInnerLoopMismtaches

## Free Energy Parameter Constants ##
R = 0.001987204258  # source: https://en.wikipedia.org/wiki/Gas_constant
T = 310.15

# Stems(source: https://rna.urmc.rochester.edu/NNDB/turner04/wc-parameters.html)
INTERMOLECULAR_INIT = 4.09  # intermolecular initiation value
STEM_SYMMETRY_PENALTY = 0.43
STEM_AU_END_PENALTY = 0.45

# Inner loops
INNER_LOOP_ASYMMETRY_PENALTY = 0.6

# other Constants
CANONICAL_BASE_PAIRS = [('A', 'U'), ('U', 'A'), ('G', 'C'),
                        ('C', 'G'), ('G', 'U'), ('U', 'G')]


class InternalLoop:
    # __init__ method for InternalLoop object
    def __init__(self, pLabel=None, label5p=None, label3p=None,  loop5p='',
                 loop3p='', loop5pSpan=(-1, -1), loop3pSpan=(-1, -1),
                 closingPairs=(('', ''), ('', '')),
                 closingPairsSpan=((-1, -1), (-1, -1)),
                 neighbor5p=('', ''), neighbor3p=('', '')):
        self._parentLabel = pLabel
        self._5pLabel = label5p
        self._3pLabel = label3p
        self._5pLoop = loop5p
        self._3pLoop = loop3p
        self._loopsLen = (len(loop5p), len(loop3p))
        self._span5p = loop5pSpan
        self._span3p = loop3pSpan
        self._closingPairs = closingPairs
        self._closingPairsSpan = closingPairsSpan
        self._neighbor5p = neighbor5p
        self._neighbor3p = neighbor5p
        self._strict = True  # used for to control energy function

    ###
    # Internal Methods
    ###

    # defines the string representation of the object
    def __str__(self):
        return f'Inner Loop: {self._parentLabel}'

    # define len function operation for InnerLoop objects
    def __len__(self):
        return self._loopsLen

    # function to update loop lengths upon change
    def _updateLoopLen(self):
        self._loopsLen = (len(self._5pLoop), len(self._3pLoop))

    # internal method to set the 5' and 3' neighbors for a InternalLoop
    def _addNeighbors(self, neighbor5p, neighbor3p):
        self._neighbor5p = neighbor5p
        self._neighbor3p = neighbor3p

    # Function checks that the internal loop has the same 5' closing pair
    # structures
    def _same5pNeighbors(self):
        return (self._neighbor5p[0] == self._neighbor5p[1])

    # Function checks that the internal loop has the same 3' closing pair
    # structures
    def _same3pNeighbors(self):
        return (self._neighbor3p[0] == self._neighbor3p[1])

    ###
    # User Accessible Methods
    ###

    '''
    Function: InternalLoop.label()
    Description: Function returns the parent label for the inner loop. Also
                allows user to set new label
    Parameters:
            (newLabel=None) -- str -- new label to identify the InternalLoop
            object
    Return Value:
            str - the current label for the InternaLoop object
    '''

    def label(self, newLabel=None):
        if newLabel:
            self._parentLabel = newLabel
        else:
            return self._parentLabel

    '''
    Function: InternalLoop.loops()
    Description: Function returns a tuple containing the sequences for the two
                inner loop subcomponents. Also allows user to define new loop
                sequences.
    Parameters:
            (loop5p=None) -- str -- new sequence to define the 5' portion of
            the internal loop
            (loop3p=None) -- str -- new sequence to define the 3' portion of
            the internal loop
    Return Value:
            (str, str) - tuple containing the 5' and 3' portions of the
            internal loop sequence
    '''

    def loops(self, loop5p=None, loop3p=None):
        if(loop5p and loop3p):
            self._5pLoop = loop5p
            self._3pLoop = loop3p
            self._updateLoopLen()
        else:
            return (self._5pLoop, self._3pLoop)

    '''
    Function: InternalLoop.loop5p()
    Description: Function that returns the 5' portion of the inner loop.
                Also allows user to define 5' portion of the loop
    Parameters:
            (loop=None) -- str -- new nucleotide sequence to define the 5'
            portion of the Internal Loop
    Return Value:
            str - the current sequence that defines the 5' portion of the
            InternalLoop
    '''

    def loop5p(self, loop=None):
        if(loop):
            self._5pLoop = loop
            self._updateLoopLen()
        else:
            return self._5pLoop

    '''
    Function: InternalLoop.loop3p()
    Description: Function that returns the 3' portion of the inner loop.
                Also allows user to define 3' portion of the loop
    Parameters:
            (loop=None) -- str -- new nucleotide sequence to define the 3'
             portion of the Internal Loop
    Return Value:
            str - the current sequence that defines the 3' portion of the
            InternalLoop
    '''

    def loop3p(self, loop=None):
        if(loop):
            self._3pLoop = loop
            self._updateLoopLen()
        else:
            return self._3pLoop

    '''
    Function: InternalLoop.loopsLen()
    Description: Function returns a tuple containing the the integer value
                lengths of the two inner loop components
    Parameters: None
    Return Value:
            (int, int) - tuple containing the integer value sequence
            lengths of the 5' and 3' portions of the Internal Loop
    '''

    def loopsLen(self):
        return self._loopsLen

    '''
    Function: InternalLoop.span()
    Description: Function returns a tuple that contains two tuples containing
                the integer start and stop positions of the 5' and 3' inner
                loop components
    Parameters: None
    Return Value:
            ((int, int), (int, int)) - tuple containing 2 tuples that each
            define the start and stop indices for the 5' and 3' portions of the
            Internal Loop
    '''

    def span(self):
        return (self._span5p, self._span3p)

    '''
    Function: InternalLoop.closingPairs()
    Description: Function returns a tuple that contains two tuples containing
                 the closing base pairs of the inner loop components
    Parameters: None
    Return Value:
            ((str, str), (str, str)) - tuple containg 3 tuples that define the
            closing nucleotide base pairs for the 5' and 3' end of the
            Internal Loop
    '''

    def closingPairs(self):
        return self._closingPairs

    '''
    Function: InternalLoop.closingPairsSpan()
    Description: Function returns a tuple that contains two tuples containing
                 the index locations of the closing base pairs of the inner
                 loop components
    Parameters: None
    Return Value:
            ((int, int), (int, int)) - uple containg 3 tuples that define the
            index locations of the closing nucleotide base pairs for the 5'
            and 3' end of the Internal Loop
    '''

    def closingPairsSpan(self):
        return self._closingPairsSpan

    '''
    Function: InternalLoop.neighbors()
    Description: function to get the StructureComponents directly adjacent
                to the InternalLoop
    Parameters: None
    Return Value:
            ((str, str), (str, str)) - tuple containg 2 tuples that define
            the neighboring structures to the 5' and 3' ends of the
            Internal Loop
    '''

    def neighbors(self):
        return (self._neighbor5p, self._neighbor3p)

    '''
    Function: InternalLoop.canonical()
    Description: Function to check if valid parameters are available to
                calculate inner loop energy
    Parameters: None
    Return Value:
            bool - returns True if there is a complete set of parameters
            for calculating the energy of the internal loop
    '''

    def canonical(self):
        # Check if energy value is present for 1x1 loop
        if len(self._5pLoop) == 1 and len(self._3pLoop) == 1:
            if (self._closingPairs[0], self._closingPairs[1], self._5pLoop, self._3pLoop) in InnerLoop_1x1_Energies:
                return True
            return False
        # check if energy value is present for 1x2 loop
        elif len(self._5pLoop) == 1 and len(self._3pLoop) == 2:
            if (self._closingPairs[0], self._closingPairs[1], self._5pLoop, self._3pLoop[1], self._3pLoop[0]) in InnerLoop_1x2_Energies:
                return True
            return False
        # check if energy value is present for 2x1 loop
        elif len(self._5pLoop) == 2 and len(self._3pLoop) == 1:
            if ((self._closingPairs[1][1], self._closingPairs[1][0]), (self._closingPairs[0][1], self._closingPairs[0][0]), self._3pLoop, self._5pLoop[1], self._5pLoop[0]) in InnerLoop_1x2_Energies:
                return True
            return False
        # Check if energy value is present for 2x2 loop
        elif len(self._5pLoop) == 2 and len(self._3pLoop) == 2:
            # convert loop sequences to proper format for dictionary
            loops = list(zip(list(self._5pLoop), list(self._3pLoop[::-1])))
            if (self._closingPairs[0], self._closingPairs[1], loops[0], loops[1]) in InnerLoop_2x2_Energies:
                return True
            return False
        # Check for valid parameters needed to calculate energy for loops
        # of other lengths
        else:
            if(self._getInnerLoopMismtachEnergy() != None):
                return True
            return False

    '''
    Function Name: _getInnerLoopInitEnergy(self)
    Description: Internal method to get the initiation energy parameter for
                Inner Loop energy function
    Parameters: None
    Return Type: float
    '''

    def _getInnerLoopInitEnergy(self):
        # get total length of inner loop for initiation parameter calculation
        loopLength = len(self._5pLoop) + len(self._3pLoop)
        if loopLength in InternalLoopInit:  # try to get initiation energy
            return float(InternalLoopInit[loopLength])
        else:  # otherwise calculate value
            return InternalLoopInit[6] + (1.08 * np.log(float(loopLength)/6.0))

    '''
    Function Name: _getInnerLoopAsymmetryEnergy(self)
    Description: Internal method to get asymmetry penalty for inner loop energy
     function
    Parameters: None
    Return Type: float
    '''

    def _getInnerLoopAsymmetryEnergy(self):
        return abs(len(self._5pLoop) - len(self._3pLoop)) * INNER_LOOP_ASYMMETRY_PENALTY

    '''
    Function Name: _getInnerLoopClosingPenalty(self)
    Description: Internal method to get the AU/GU Closing penalty for
                InnerLoop energy function
    Parameters: None
    Return Type: float
    '''

    def _getInnerLoopClosingPenalty(self):
        closingPenalty = 0
        # closing pairs that result in end penalty
        endPenaltyPairs = [('A', 'U'), ('G', 'U'), ('U', 'A'), ('U', 'G')]
        # get the closing pairs for the inner loop
        closingPair5p, closingPair3p = self.closingPairs()
        # check for penalty condition in 5' closing pair
        if closingPair5p in endPenaltyPairs:
            closingPenalty += 0.7
        # check for penalty in 3' closing pair
        if closingPair3p in endPenaltyPairs:
            closingPenalty += 0.7

        return float(closingPenalty)

    '''
    Function Name: _getInnerLoopMismatchEnergy_3x2(self)
    Description: Internal method to get the mismatch energy for a 3x2 InnerLoop
    Parameters: None
    Return Type: float
    '''

    def _getInnerLoopMismatchEnergy_3x2(self):
        loop1, loop2 = self.loops()
        mismatch5p = (loop2[0], loop1[-1])
        mismatch3p = (loop1[0], loop2[-1])

        mismatchEnergy_3x2 = 0
        # check for mismatch condition between 5' closing pair and first mismatch
        if ((self._closingPairs[1][1], self._closingPairs[1][0]), mismatch5p) in InnerLoopMismatches_2x3:
            mismatchEnergy_3x2 += InnerLoopMismatches_2x3[(
                (self._closingPairs[1][1], self._closingPairs[1][0]), mismatch5p)]
        else:
            if (self._strict):
                return None

        # check for mismatch condition between 3'closing pair and mismatch 2
        if ((self._closingPairs[0][1], self._closingPairs[0][0]), mismatch3p) in InnerLoopMismatches_2x3:
            mismatchEnergy_3x2 += InnerLoopMismatches_2x3[(
                (self._closingPairs[0][1], self._closingPairs[0][0]), mismatch3p)]
        else:
            if (self._strict):
                return None

        return float(mismatchEnergy_3x2)

    '''
    Function Name: _getInnerLoopMismatchEnergy_2x3(self)
    Description: Internal method to get the mismatch energy for a 2x3 InnerLoop
    Parameters: None
    Return Type: float
    '''

    def _getInnerLoopMismatchEnergy_2x3(self):
        loop1, loop2 = self.loops()  # get both loops
        mismatch5p = (loop1[0], loop2[-1])  # get 1st mismatch
        mismatch3p = (loop2[0], loop1[-1])  # get 2nd mismatch

        mismatchEnergy_2x3 = 0
        # check for mismatch condition between 5' closing pair and first mismatch
        if (self._closingPairs[0], mismatch5p) in InnerLoopMismatches_2x3:
            mismatchEnergy_2x3 += InnerLoopMismatches_2x3[(
                self._closingPairs[0], mismatch5p)]
        else:
            if (self._strict):
                return None

        # check for mismatch condition between 3'closing pair and mismatch 2
        if ((self._closingPairs[1][1], self._closingPairs[1][0]), mismatch3p) in InnerLoopMismatches_2x3:
            mismatchEnergy_2x3 += InnerLoopMismatches_2x3[(
                (self._closingPairs[1][1], self._closingPairs[1][0]), mismatch3p)]
        else:
            if (self._strict):
                return None

        return float(mismatchEnergy_2x3)

    '''
    Function Name: _getInnerLoopMismatchEnergy_Other(self)
    Description: Internal method to get the inner loop mismtach energy
                for other inner loops
    Parameters: None
    Return Type: float
    '''

    def _getInnerLoopMismatchEnergy_Other(self):
        loop1, loop2 = self.loops()  # get both loops
        mismatch5p = (loop1[0], loop2[-1])  # get 1st mismatch
        mismatch3p = (loop1[-1], loop2[0])  # get 2nd mismatch

        mismatchEnergy_Other = 0
        # check for mismatch 1 for condition
        if mismatch5p in OtherInnerLoopMismtaches:
            mismatchEnergy_Other += OtherInnerLoopMismtaches[mismatch5p]
        elif (self._strict):
            return None

        # check mismatch 2 for condition
        if mismatch3p in OtherInnerLoopMismtaches:
            mismatchEnergy_Other += OtherInnerLoopMismtaches[mismatch3p]
        elif (self._strict):
            return None

        return float(mismatchEnergy_Other)

    '''
    Function Name: _getInnerLoopMismtachEnergy(self)
    Description: Internal method to get the mismatch energy for an inner loop
    Parameters: None
    Return Type: float
    '''

    def _getInnerLoopMismtachEnergy(self):
        # 1 x (n-1) Inner Loops
        loopLength = len(self._5pLoop) + len(self._3pLoop)
        if (len(self._5pLoop) == 1 and len(self._3pLoop) == loopLength-1) or (len(self._5pLoop) == loopLength-1 and len(self._3pLoop) == 1):
            return 0.0  # Mismatch energy is 0, so we dont need to do anything

        # 2x3 Inner Loop mismatches
        elif (len(self._5pLoop) == 2 and len(self._3pLoop) == 3):
            return self._getInnerLoopMismatchEnergy_2x3()

        # 3x2 inner loop mismatches
        elif (len(self._5pLoop) == 3 and len(self._3pLoop) == 2):
            return self._getInnerLoopMismatchEnergy_3x2()

        # other inner loops
        else:
            return self._getInnerLoopMismatchEnergy_Other()

    '''
    Function Name: _calcEnergy(self)
    Description: Internal method  to calculate the energy for inner loops
                whose energies are not stored in the imported dictionaries
    Parameters: None
    Return Type: float
    '''

    def _calcEnergy(self):
        # get InnerLoop initiation parameter
        ilInit = self._getInnerLoopInitEnergy()
        if(ilInit is None):  # check that parameter is present
            return None

        # asymmetry penalty
        asym = self._getInnerLoopAsymmetryEnergy()
        if(asym is None):  # check that parameter is present
            return None

        # AU / GU Closure penalty
        closingPenalty = self._getInnerLoopClosingPenalty()
        if(closingPenalty is None):  # check that parameter is present
            return None

        # get mismtach energy
        mismatchEnergy = self._getInnerLoopMismtachEnergy()
        if(mismatchEnergy is None):  # check that parameter is present
            return None

        # sum energy components and return
        return ilInit + asym + closingPenalty + mismatchEnergy

    '''
    Function Name: energy(self)
    Description: Function to get the free energy for the inner loop object
    Parameters: None
    Return Type: float
    '''

    def energy(self, strict=True):
        # set mode for energy calculations
        self._strict = strict

        # check for 1x1 - value taken from imported dicitionary
        if len(self._5pLoop) == 1 and len(self._3pLoop) == 1:
            # check if key in dictionary
            if (self._closingPairs[0], self._closingPairs[1], self._5pLoop, self._3pLoop) in InnerLoop_1x1_Energies:
                loopEnergy = InnerLoop_1x1_Energies[(
                    self._closingPairs[0], self._closingPairs[1], self._5pLoop, self._3pLoop)]
                return loopEnergy
            else:  # otherwise calculate energy
                if(self._strict):
                    return None
                else:
                    return self._calcEnergy()

        # check for 1x2 - value taken from imported dicitionary
        elif len(self._5pLoop) == 1 and len(self._3pLoop) == 2:
            # check if key in dictionary
            if (self._closingPairs[0], self._closingPairs[1], self._5pLoop, self._3pLoop[1], self._3pLoop[0]) in InnerLoop_1x2_Energies:
                loopEnergy = InnerLoop_1x2_Energies[(
                    self._closingPairs[0], self._closingPairs[1], self._5pLoop, self._3pLoop[1], self._3pLoop[0])]
                return loopEnergy
            else:  # otherwise calculate energy
                if(self._strict):
                    return None
                else:
                    return self._calcEnergy()

        # check for 2x1 case - value taken from dicitonary
        elif len(self._5pLoop) == 2 and len(self._3pLoop) == 1:
            if ((self._closingPairs[1][1], self._closingPairs[1][0]), (self._closingPairs[0][1], self._closingPairs[0][0]), self._3pLoop, self._5pLoop[1], self._5pLoop[0]) in InnerLoop_1x2_Energies:  # check if key in dictionary
                loopEnergy = InnerLoop_1x2_Energies[((self._closingPairs[1][1], self._closingPairs[1][0]), (
                    self._closingPairs[0][1], self._closingPairs[0][0]), self._3pLoop, self._5pLoop[1], self._5pLoop[0])]
                return loopEnergy
            else:  # otherwise calculate energy
                if(self._strict):
                    return None
                else:
                    return self._calcEnergy()

        # check for 2x2 - value taken from imported dicitionary
        elif len(self._5pLoop) == 2 and len(self._3pLoop) == 2:
            # convert loop sequences to proper format for dictionary
            loops = list(zip(list(self._5pLoop), list(self._3pLoop[::-1])))
            # check if key in dictionary
            if (self._closingPairs[0], self._closingPairs[1], loops[0], loops[1]) in InnerLoop_2x2_Energies:
                loopEnergy = InnerLoop_2x2_Energies[(
                    self._closingPairs[0], self._closingPairs[1], loops[0], loops[1])]
                return loopEnergy
            else:  # otherwise calculate energy
                if(self._strict):
                    return None
                else:
                    return self._calcEnergy()

        # Other cases need to be calculated
        else:
            return self._calcEnergy()
