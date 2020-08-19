'''
Filename: Structure.py
Author: Michael Hathaway

Description: python module that defines the Structure Object.
The Structure Object provides a user friendly mechanism for working with
RNA structure type files in the python programming language.
'''


'''
## About the structure object ##
The Structure object is a python object-oriented representation of the
information contained within an RNA Structure Type file. The object provides a
mechanism to easily access and work with the data in the python programming
language. In addition it includes functionality for calculating the energy
associated with certain RNA secondary structures within the RNA molecule.
'''


class Structure:
    # __init__() method for the Structure object
    def __init__(self, filename=None):
        # RNA Molecule basic info
        # all values are stored as strings
        self._name = None
        self._length = None
        self._pageNum = None

        # structural representations of RNA molecule
        # all values are stored as strings
        self._sequence = None
        self._DBN = None
        self._structureArray = None
        self._varna = None

        '''
        secondary structure information
        secondary structure information for each RNA molecule is stored as a
        dictionary where the key is the label for the secondary structure and
        the value is an object that contains the data for the secondary
        structure with the given label.

        For information on how each secondary structure class was implemented
        and how to access secondary structure data, see the
        StructureTypeComponents section below.
        '''
        self._stems = {}
        self._hairpins = {}
        self._bulges = {}
        self._internalLoops = {}
        self._multiLoops = {}
        self._externalLoops = {}
        self._pk = {}
        self._ncbp = {}
        self._ends = {}

        '''
        Component Array
        The component array is a numpy array of the same length as the molecule
        where each index contains the label for the secondary structure tha
        that index is a part of. the component array is initialized as None.
        When the length of the molecule is parsed from the .st file, a numpy
        array of that length is generated
        '''
        self._componentArray = None

    # define string representation of the molecule
    def __str__(self):
        return f'RNA: {self._name}'

    # define len function for Structure object
    def __len__(self):
        return self._length


####################################################
###### Load File and Associated Functions ##########
####################################################

    '''
    Function: resetStructure()
    Description: Function resets all Structure member variables to None.
    Function is called when a file cannot be completley parsed.
    Parameters:
            None
    Return Value:
            None
    '''

    def resetStructure(self):
        # reset basic molecule information
        self._name = None
        self._length = None
        self._pageNum = None

        # reset all structure representations of the molecule
        self._sequence = None
        self._DBN = None
        self._structureArray = None
        self._varna = None

        # reset all StructureComponent dictionaries
        self._stems.clear()
        self._hairpins.clear()
        self._bulges.clear()
        self._internalLoops.clear()
        self._multiLoops.clear()
        self._externalLoops.clear()
        self._pk.clear()
        self._ncbp.clear()
        self._ends.clear()

        # reset component array
        self._componentArray = None


##############################################
###### Add StructureComponent Neighbors ######
##############################################

    '''
    Function Name: _addStructureComponentNeighbors()
    Description: Function fills in the neighboring structure information for each of the StructureComponent objects contained in the Structure object
    Parameters:
            None
    Return Value:
            None
    '''

    def _addStructureComponentNeighbors(self):
        # get all labels for StructureComponents in Structure
        allStructureComponentLabels = self.features()

        for feature in allStructureComponentLabels:  # iterate through all labels
            try:
                # get the StructureComponent object for a given label
                structureComponent = self.component(feature)
                # get neighbors of the feature
                neighbors = self.neighbors(feature)
                # add neighbors to StructureComponent object
                structureComponent._addNeighbors(neighbors[0], neighbors[1])
            except:  # will skip NCBPs and MultiLoops
                continue


#############################
###### COMPONENT ARRAY ######
#############################

    '''
    Function Name: _addStemToComponentArray(stem)
    Description: Internal method used in _loadFile() that adds a given stem to the component array
    Parameters:
            (stem) - Stem object - stem to be added to the component array
    Return Type:
            None
    '''

    def _addStemToComponentArray(self, stem):
        for i in range(stem.sequence5pSpan()[0]-1, stem.sequence5pSpan()[1]):
            self._componentArray[i] = stem.label()

        for i in range(stem.sequence3pSpan()[0]-1, stem.sequence3pSpan()[1]):
            self._componentArray[i] = stem.label()

    '''
    Function Name: _addBulgeToComponentArray(bulge)
    Description:  Internal method used in _loadFile() that adds a given bulge to the component array
    Parameters:
            (bulge) - Bulge object - bulge to be added to the component array
    Return Type:
            None
    '''

    def _addBulgeToComponentArray(self, bulge):
        for i in range(bulge.span()[0]-1, bulge.span()[1]):
            self._componentArray[i] = bulge.label()

    '''
    Function Name: _addHairpinToComponentArray(hairpin)
    Description: Internal method used in _loadFile() that adds a hairpin to the component array
    Parameters:
            (hairpin) - Hairpin object - hairpin to be added to the component array
    Return Type:
            None
    '''

    def _addHairpinToComponentArray(self, hairpin):
        for i in range(hairpin.span()[0]-1, hairpin.span()[1]):
            self._componentArray[i] = hairpin.label()

    '''
    Function Name: _addEndToComponentArray(end)
    Description: Internal method used in _loadFile() that adds an end to the component array
    Parameters:
            (end) - End object - end to be added to the component array
    Return Type:
            None
    '''

    def _addEndToComponentArray(self, end):
        for i in range(end.span()[0]-1, end.span()[1],):
            self._componentArray[i] = end.label()

    '''
    Function Name: _addInternalLoopToComponentArray(InternalLoop)
    Description: Internal method used in _loadFile() that adds an inner loop to the component array
    Parameters:
            (internalLoop) - InternalLoop object - inner loop to be added to the component array
    Return Type:
            None
    '''

    def _addInternalLoopToComponentArray(self, internalLoop):
        for pair in internalLoop.span():
            for i in range(pair[0]-1, pair[1]):
                self._componentArray[i] = internalLoop.label()

    '''
    Function Name: _addExternalLoopToComponentArray(el)
    Description: Internal method used in _loadFile() that adds an external loop to the component array
    Parameters:
            (el) - ExternalLoop object - External loop to be added to the component array
    Return Type:
            None
    '''

    def _addExternalLoopToComponentArray(self, el):
        for i in range(el.span()[0]-1, el.span()[1]):
            self._componentArray[i] = el.label()

    '''
    Function Name: _addMultiLoopToComponentArray(multiloop)
    Description: Internal method used in _loadFile() that adds an multiloop to the component array
    Parameters:
            (multiloop) - MultiLoop object - multiloop to be added to the component array
    Return Type:
             None
    '''

    def _addMultiLoopToComponentArray(self, multiloop):
        for subunit in multiloop.subunitLabels():  # iterate through subunit labels
            span = multiloop._spans[subunit]  # get span for particular subunit
            for i in range(span[0]-1, span[1]):
                self._componentArray[i] = multiloop._parentLabel

    '''
    Function Name: componentArray()
    Description: function that returns the _componentArray for the Structureobject
    Parameters:
            None
    Return Type:
            numpy array
    '''

    def componentArray(self):
        return self._componentArray


###########################
###### SEQUENCE INFO ######
###########################

    '''
    Function Name: Name()
    Description: Function returns that name of the RNA molecule represented in the .st file
    Parameters:
            None
    Return Type:
            str
    '''

    def name(self):
        return self._name

    '''
    Function Name: Length()
    Description: Function returns the length of the RNA sequence represented in the .st file
    Parameters:
            None
    Return Type:
            int
    '''

    def length(self):
        return self._length

    '''
    Function Name: PageNum()
    Description: Function returns the page number for the RNA molecule represented in the .st file
    Parameters:
            None
    Return Type:
            int
    '''

    def pageNum(self):
        return self._pageNum


########################################
###### STRUCTURAL REPRESENTATIONS ######
########################################

    '''
    Function Name: Sequence()
    Description: Function returns the RNA sequence(A,U,G,C) for the RNA molecule represented in the .st file
    Parameters:
            None
    Return Type:
            str
    '''

    def sequence(self):
        return self._sequence

    '''
    Function Name: DotBracket()
    Description: function to get the Dot Bracket notation for the Structureobject
    Parameters:
            None
    Return Type:
            str
    '''

    def dotBracket(self):
        return self._DBN

    '''
    Function Name: StructureArray()
    Description: Function to get the structure Array for the StructureTyoe object
    Parameters:
            None
    Return Type:
            str
    '''

    def structureArray(self):
        return self._structureArray

    '''
    Function Name: VARNA()
    Description:
    Parameters:
            None
    Return Type:
             str
    '''

    def VARNA(self):
        return self._varna


###################
###### STEMS ######
###################

    '''
    Function Name: _addStemBulgeNeighborBooleans()
    Description: Function fills in boolean values for whether or not each stem object is adjacent to a bulge of length 1
    Parameters:
            None
    Return value:
            None
    '''

    def _addStemBulgeNeighborBooleans(self):
        for stem in self.stems():  # iterate through stems
            neighbor5p, neighbor3p = self.neighbors(
                stem.label(), object=True)  # get neighbors
            bool5p, bool3p = False, False

            # check if a 5' neighbor is a length=1 bulge and if so change bool5p to True
            if(neighbor5p[0] != 'EOM' and neighbor5p[0].label()[0] == 'B'):
                if (neighbor5p[0].sequenceLen() == 1):
                    bool5p = True
            elif(neighbor5p[1] != 'EOM' and neighbor5p[1].label()[0] == 'B'):
                if (neighbor5p[1].sequenceLen() == 1):
                    bool5p = True

            # check if a 3' neighbor is a length=1 bulge and if so change bool3p to True
            if(neighbor3p[0] != 'EOM' and neighbor3p[0].label()[0] == 'B'):
                if (neighbor3p[0].sequenceLen() == 1):
                    bool3p = True
            elif(neighbor3p[1] != 'EOM' and neighbor3p[1].label()[0] == 'B'):
                if (neighbor3p[1].sequenceLen() == 1):
                    bool3p = True

            stem._addAdjacentBulgeBoolean(bool5p, bool3p)

    '''
    Function Name: addStem()
    Description: Function to add a new stem to the Structureobject
    Parameters:
            (stemLabel) - str - key for stem object in self._stems dictionary
            (newStem) - Stem object - Stem object to be stored at given key in the self._stems dictionary
    Return Value:
            None
    '''

    def addStem(self, stemLabel, newStem):
        self._stems[stemLabel] = newStem

    '''
    Function Name: stemLabels()
    Description: function to access all the stem labels for the Structureobject
    Parameters:
            None
    Return Type:
            list
    '''

    def stemLabels(self):
        return list(self._stems.keys())

    '''
    Function Name: stems(label=None)
    Description: function to get all the stem objects in the Structureobject. If label is provided the stem object
    with the matching label is returned.
    Parameters:
            (label=None) - str - the label for the stem being searched for.
    Return Type:
            list or Stem object
    '''

    def stems(self, label=None):
        if label:
            return self._getStemByLabel(label)
        else:
            return list(self._stems.values())

    '''
    Function Name: numStems()
    Description: function to get the number of stems in a Structureobject
    Parameters:
             None
    Return Value:
             int
    '''

    def numStems(self):
        return len(self._stems)

    '''
    Function Name: getStemByLabel(stemLabel)
    Description: Function to get a particular Stem object based on its label
    Parameters:
            (stemLabel) - str - label for Stem to be accessed
    Return Value:
             Stem object
    '''

    def _getStemByLabel(self, stemLabel):
        try:
            stem = self._stems[stemLabel]
            return stem
        except KeyError:
            print(f'Stem label: {stemLabel} not found.')
            return None


######################
###### HAIRPINS ######
######################

    '''
    Function Name: addHairpin(label, newHairpin)
    Description: Function to add a new hairpin to the Structureobject
    Parameters:
            (label) - str - key for Hairpin object in self._hairpins dictionary
            (newHairpin) - Hairpin object - Hairpin object to be stored at the given key in the self._haripins dicitonary
    Return Type:
            None
    '''

    def addHairpin(self, label, newHairpin):
        self._hairpins[label] = newHairpin

    '''
    Function Name: hairpinLabels()
    Description: Function to access all the hairpin labels in the Structureobject
    Parameters:
            None
    Return Type:
             list
    '''

    def hairpinLabels(self):
        return list(self._hairpins.keys())

    '''
    Function Name: hairpins()
    Description: Function to get all Hairpin objects in the Structureobject. if label is provided, the Hairpin object with the
    matching label is returned.
    Parameters:
            None
    Return Type:
             list
    '''

    def hairpins(self, label=None):
        if label:
            return self._getHairpinByLabel(label)
        else:
            return list(self._hairpins.values())

    '''
    Function Name: numHairpins()
    Description: Function to get the number of hairpins in the Structureobject
    Parameters:
            None
    Return Type:
             int
    '''

    def numHairpins(self):
        return len(self._hairpins)

    '''
    Function Name: getHairpinByLabel(hairpinLabel)
    Description: Function to access a particular Hairpin Object based on its label
    Parameters:
            (hairpinLabel) - str - label for hairpin to accessed
    Return Type:
             Hairpin Object
    '''

    def _getHairpinByLabel(self, hairpinLabel):
        try:
            hairpin = self._hairpins[hairpinLabel]
            return hairpin
        except KeyError:
            print(f'Hairpin label: {hairpinLabel} not found')


####################
###### BULGES ######
####################

    '''
    Function Name: addBulge(bulgeLabel, newBulge)
    Description: Function to add a new Bulge object to the Structureobject
    Parameters:
            (bulgeLabel) - str - the key value to be used for the new Bulge object
            (newBulge) - Bulge Object - Bulge object to be stored at the given key in the self._bulges dictionary
    Return Type:
             None
    '''

    def addBulge(self, bulgeLabel, newBulge):
        self._bulges[bulgeLabel] = newBulge

    '''
    Function Name: bulgeLabels()
    Description: Function to get all the bulge labels for the Structureobject
    Parameters:
            None
    Return Type:
             list
    '''

    def bulgeLabels(self):
        return list(self._bulges.keys())

    '''
    Function Name: bulges()
    Description: Function to get all the Bulge objects for the Structureobject
    Parameters:
            None
    Return Type:
             list
    '''

    def bulges(self, label=None):
        if label:
            return self._getBulgeByLabel(label)
        else:
            return list(self._bulges.values())

    '''
    Function Name: numBulges()
    Description: Function to get the number of bulges in a given Structureobject
    Parameters:
            None
    Return Type:
             int
    '''

    def numBulges(self):
        return len(self._bulges)

    '''
    Function Name: getBulgeByLabel(bulgeLabel)
    Description: Function to access a particular Bulge object based on its label
    Parameters:
            (bulgeLabel) - str - label for Bulge object to be accessed
    Return Type:
             Bulge object
    '''

    def _getBulgeByLabel(self, bulgeLabel):
        try:
            bulge = self._bulges[bulgeLabel]
            return bulge
        except KeyError:
            print(f'Bulge label: {bulgeLabel} not found')
            return None


###########################
###### Inner Loops ########
###########################

    '''
    Function Name: addInternalLoop(parentLabel, subunitLabel, newInternalLoop)
    Description: function to add a new InternalLoop object to the Structureobject
    Parameters:
            (parentLabel) - str - parent key value for the Inner loop to be added
            (newInternalLoop) - InternalLoop object - InternalLoop object to be stored at the the given key in the self._internalLoops dictionary
    Return Type:
            None
    '''

    def addInternalLoop(self, parentLabel, newInternalLoop):
        self._internalLoops[parentLabel] = newInternalLoop

    '''
    Function Name: InternalLoopLabels()
    Description: Function to return a list of all the inner loop labels in the Structure object
    Parameters:
            None
    Return Type:
             list
    '''

    def internalLoopLabels(self):
        return list(self._internalLoops.keys())

    '''
    Function Name: InternalLoops()
    Description: Function to return a list of the InternalLoop objects in the Structure object
    Parameters:
            None
    Return Type:
            list
    '''

    def internalLoops(self, label=None):
        if label:
            return self._getInternalLoopByLabel(label)
        else:
            return list(self._internalLoops.values())

    '''
    Function Name: numInternalLoops()
    Description: function to get the number of inner loops in a Structure object
    Parameters:
            None
    Return Type:
             int
    '''

    def numInternalLoops(self):
        return len(self._internalLoops)

    '''
    Function Name: getInternalLoopByLabel(label)
    Description: returns InternalLoop object stored at the given key value
    Parameters:
            (label) - str - the key value for the inner loop to be accessed
    Return Type:
             InternalLoop object
    '''

    def _getInternalLoopByLabel(self, label):
        try:
            internalLoop = self._internalLoops[label]
            return internalLoop
        except KeyError:
            print(f'Internal Loop: {label} not found.')
            return None

    '''
    Function Name: getInternalLoopSubunitByLabel(parentLabel, subunitLabel)
    Description:
    Parameters:
             (parentLabel) - str - parent label for inner loop object
            (subunitLabel) - str - subunit label for the inner loop object
    Return Type:
            dictionary with the following subunit information

                {
                    'label' : ,
                    'Sequence' : ,
                    'span' :
                }
    '''

    def getInternalLoopSubunit(self, parentLabel, subunitLabel):
        internalLoop = None
        try:
            internalLoop = self._internalLoops[parentLabel]
        except KeyError:
            print(f'Inner Loop: {parentLabel} not found.')

        if subunitLabel == '1':
            subunit = {
                'label': f'{internalLoop._parentLabel}.1',
                'Sequence': internalLoop._5pSequence,
                'span': internalLoop._span5p,
            }
            return subunit
        elif subunitLabel == '2':
            subunit = {
                'label': f'{internalLoop._parentLabel}.2',
                'Sequence': internalLoop._3pSequence,
                'span': internalLoop._span5p,
            }
            return subunit
        else:
            return None


###########################
###### MultiLoops #########
###########################

    '''
    Function Name: addMultiLoop(parentLabel, subunitLabel, newMultiLoop)
    Description: function to add a new MultiLoop object to the Structure object
    Parameters:
            (parentLabel) - str - parent multiloop label for the multiloop to be added
            (subunitLabel) - str - subunit label for the multiloop to be added
            (newMultiLoop) - MultiLoop object - MultiLoop object to be added at the given key values
    Return Type:
             None
    '''

    def addMultiLoop(self, parentLabel, newMultiLoop):
        self._multiLoops[parentLabel] = newMultiLoop

    '''
    Function Name: numMultiLoops()
    Description: Function to get the number of multiloops in a StructureTyoe object
    Parameters:
            None
    Return Type:
             int
    '''

    def numMultiLoops(self):
        return len(self._multiLoops)

    '''
    Function Name: _getMultiLoopByLabel(self, label)
    Description: Function returns the MultiLoop object identified by a given parent label.
    Parameters:
            (label) -- str -- parent label of the MultiLoop object being accessed
    Return Value:
            MultiLoop object
    '''

    def _getMultiLoopByLabel(self, label):
        try:
            multiloop = self._multiLoops[label]
            return multiloop
        except KeyError:
            print(f'MultiLoop {label} not found.')
            return None

    '''
    Function Name: multiLoops(label=None)
    Description: Function returns a list of all the multiloops stored in a Structure object. If a particular label is provided, only the MultiLoop object with that label will be returned
    Parameters:
            (label=None) -- str -- parent label of the MultiLoop object being accessed. Default value is None
    Return Value:
            if label is provided: MultiLoop object
            if label=None: list of MultiLoop objects
    '''

    def multiLoops(self, label=None):
        if(label):
            return self._getMultiLoopByLabel(label)
        else:
            return list(self._multiLoops.values())


###############################
###### External Loops #########
###############################

    '''
    Function Name: addExternalLoop(elLabel, newEL)
    Description: function to add a new External Loop object to the Structure object
    Parameters:
            (elLabel) - str - the key value to be used for the new ExternalLoop object
            (newEL) - ExternalLoop Object - External loop to be stored at given key value
    Return Type:
             None
    '''

    def addExternalLoop(self, elLabel, newEL):
        self._externalLoops[elLabel] = newEL

    '''
    Function Name: externalLoopLabels()
    Description: Function to return a list of the external loop labels for a Structure object
    Parameters:
            none
    Return Type:
             list
    '''

    def externalLoopLabels(self):
        return list(self._externalLoops.keys())

    '''
    Function Name: externalLoops()
    Description: function to return a list of all the ExternalLoop Objects in a Structure object
    Parameters:
            None
    Return Type:
             list
    '''

    def externalLoops(self, label=None):
        if label:
            return self._getExternalLoopByLabel(label)
        else:
            return list(self._externalLoops.values())

    '''
    Function Name: numExternalLoops()
    Description: function to return the number of external loops in a Structure object
    Parameters:
            None
    Return Type:
             int
    '''

    def numExternalLoops(self):
        return len(self._externalLoops)

    '''
    Function Name: getExternalLoopByLabel(elLabel)
    Description: Function to access a particular ExternalLoop Object based on its label
    Parameters:
            (elLabel) - str - label for the ExternalLoop to be accessed
    Return Type:
             ExternalLoop object
    '''

    def _getExternalLoopByLabel(self, elLabel):
        try:
            el = self._externalLoops[elLabel]
            return el
        except KeyError:
            print(f'External Loop: {elLabel} not found.')
            return None


####################
###### NCBP ########
####################

    '''
    Function Name: addNCBP(ncbpLabel, newNCBP)
    Description: Function to add a new NCBP to the Structure object
    Parameters:
            (ncbpLabel) - str - label for the new NCBP
            (newNCBP) - NCBP Object - NCBP object to be added
    Return Type:
             None
    '''

    def addNCBP(self, ncbpLabel, newNCBP):
        self._ncbp[ncbpLabel] = newNCBP

    '''
    Function Name: ncbpLabels()
    Description: function to return a list of all the ncbp labels in a given Structure object
    Parameters:
            None
    Return Type:
             list
    '''

    def ncbpLabels(self):
        return list(self._ncbp.keys())

    '''
    Function Name: NCBPs()
    Description: function to return a list of all the NCBPs in a Structure object
    Parameters:
            None
    Return Type:
             list
    '''

    def NCBPs(self, label=None):
        if label:
            return self._getNCBPByLabel(label)
        else:
            return list(self._ncbp.values())

    '''
    Function Name: numNCBPs()
    Description: Function to get the number of NCBP's in a given Structure object
    Parameters:
            None
    Return Type:
             int
    '''

    def numNCBPs(self):
        return len(self._ncbp)

    '''
    Function Name: getNCBPByLabel(ncbpLabel)
    Description: Function to get a particular NCBP object based on its label
    Parameters:
            (ncbpLabel) - str - label for NCBP object to be accessed
    Return Type:
             NCBP object
    '''

    def _getNCBPByLabel(self, ncbpLabel):
        try:
            ncbp = self._ncbp[ncbpLabel]
            return ncbp
        except KeyError:
            print(f'NCBP label: {ncbpLabel} not found.')
            return None


####################
###### ENDS ########
####################

    '''
    Function Name: addEnd(endLabel, newEnd)
    Description: Function to add a new End to the Structure object
    Parameters:
            (endLabel) - str - label for new End object
            (newEnd) - End Object - new End object to be added
    Return Type:
             None
    '''

    def addEnd(self, endLabel, newEnd):
        self._ends[endLabel] = newEnd

    '''
    Function Name: endLabels()
    Description: Function to return a list of all the end labels for the Structure object
    Parameters:
            None
    Return Type:
             list
    '''

    def endLabels(self):
        return list(self._ends.keys())

    '''
    Function Name: ends()
    Description: function to return a list of all the End objects for the Structure object
    Parameters:
            (label=None) - str - label for the end being accessed
    Return Type:
            End object or list
    '''

    def ends(self, label=None):
        if label:
            return self._getEndByLabel(label)
        else:
            return list(self._ends.values())

    '''
    Function Name: numEnds()
    Description: function to get the number of ends in Structure object
    Parameters:
            None
    Return Type:
             int
    '''

    def numEnds(self):
        return len(self._ends)

    '''
    Function Name: getEndByLabel(endLabel)
    Description: function to access a particular End Object based on its label
    Parameters:
            (endLabel) - str - the label for the End object to be accessed
    Return Type:
             End Object
    '''

    def _getEndByLabel(self, endLabel):
        try:
            end = self._ends[endLabel]
            return end
        except KeyError:
            print(f'End: {endLabel} not found.')
            return None


################################
######## OTHER FUNCTIONs #######
################################

    '''
    Function Name: component(label, subLabel=None)
    Description: More general form of get<secondarystructure>ByLabel(). Allows you to access a given
    structure type component based on its label. Useful when using the componentArray because you may not know which
    type of structure type component you will be accessing if you are(for example) looping through the entire array
    Parameters:
            (label) - str - label of the feature to be accessed
    Return Type:
            StructureComponent object
    '''

    def component(self, label):
        if label[0] == 'S':
            return self._getStemByLabel(label)
        elif label[0] == 'H':
            return self._getHairpinByLabel(label)
        elif label[0] == 'B':
            return self._getBulgeByLabel(label)
        elif label[0] == 'X':
            return self._getExternalLoopByLabel(label)
        elif label[0] == 'E':
            return self._getEndByLabel(label)
        elif label[0] == 'N':
            return self._getNCBPByLabel(label)
        elif label[0] == 'I':
            return self._getInternalLoopByLabel(label)
        elif label[0] == 'M':
            return self._getMultiLoopByLabel(label)
        else:
            # if label is not handled by any of these blocks
            print(f'Label: {label} not found in Structure object.')
            return None

    '''
    Function Name: neighbors(label, object=False)
    Description: Function to get the secondary structures adjacent to the feature of interest.
    Parameters:
            (label) - str - label for the feature of interest
            (object) - bool - optional argument that causes the function to return the actual StructureTypeComponent objects instead of just the object label
    Return Type:
            Returns a tuple containing the labels for the adjacent features in order of 5' to 3' locations
    '''

    def neighbors(self, label, object=False):
        adjacentFeatures = []  # list to store the adjacent RNA features

        if label in self._componentArray:  # check if the feature is valid
            # get index locations of the feature
            span = self.component(label).span()
            # tuple only containes integer index locations(example: bulge location)
            if all(type(i) is int for i in span):
                try:  # try/except block will handle ends which only have one neighbor and one out of range index
                    neighbor5p = (self._componentArray[span[0]-2] if not object else self.component(
                        self._componentArray[span[0]-2]))
                except:
                    neighbor5p = 'EOM'  # 'End of Molecule'

                try:  # try/except block will handle ends which only have one neighbor and one out of range index
                    neighbor3p = (self._componentArray[span[1]] if not object else self.component(
                        self._componentArray[span[1]]))
                except:
                    neighbor3p = 'EOM'

                adjacentFeatures = (neighbor5p, neighbor3p)

            # tuple containes other tuples within it(example: InternalLoop locations)
            else:
                try:  # try/except block will handle ends which only have one neighbor and one out of range index
                    seq1_neighbor5p = (
                        self._componentArray[span[0][0]-2] if not object else self.component(self._componentArray[span[0][0]-2]))
                except:
                    seq1_neighbor5p = 'EOM'

                try:  # try/except block will handle ends which only have one neighbor and one out of range index
                    seq1_neighbor3p = (self._componentArray[span[0][1]] if not object else self.component(
                        self._componentArray[span[0][1]]))
                except:
                    seq1_neighbor3p = 'EOM'

                try:  # try/except block will handle ends which only have one neighbor and one out of range index
                    seq2_neighbor5p = (
                        self._componentArray[span[1][0]-2] if not object else self.component(self._componentArray[span[1][0]-2]))
                except:
                    seq2_neighbor5p = 'EOM'

                try:  # try/except block will handle ends which only have one neighbor and one out of range index
                    seq2_neighbor3p = (self._componentArray[span[1][1]] if not object else self.component(
                        self._componentArray[span[1][1]]))
                except:
                    seq2_neighbor3p = 'EOM'

                adjacentFeatures = (
                    (seq1_neighbor5p, seq2_neighbor3p), (seq2_neighbor5p, seq1_neighbor3p))

            return adjacentFeatures
        else:  # otherwise print error and return None
            return None

    """
    Function: features()
    Description: Function to return a list of all the StructureComponent labels in a bpRNAStructure object
    Parameters:
            None
    Return Value:
            list of strings, each string is the label for a StructureComponent in the bpRNAStructure
            *Currently only return the StructureComponents: stems, bulges, hairpins, and inner loops
    """

    def features(self):
        featureLabels = []

        featureLabels.extend(self.stemLabels())
        featureLabels.extend(self.bulgeLabels())
        featureLabels.extend(self.hairpinLabels())
        featureLabels.extend(self.internalLoopLabels())

        return featureLabels
