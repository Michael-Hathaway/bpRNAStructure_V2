import numpy as np
import re

# Structure Type Imports
from bpRNAStructure.Structure import Structure

from bpRNAStructure.StructureComponents.Hairpin import Hairpin
from bpRNAStructure.StructureComponents.Stem import Stem
from bpRNAStructure.StructureComponents.Bulge import Bulge
from bpRNAStructure.StructureComponents.InternalLoop import InternalLoop
from bpRNAStructure.StructureComponents.ExternalLoop import ExternalLoop
from bpRNAStructure.StructureComponents.Multiloop import Multiloop
from bpRNAStructure.StructureComponents.End import End
from bpRNAStructure.StructureComponents.NCBP import NCBP

from bpRNAStructure.StructureExecptions import StructureFileNotProvided


class StructureCreator():
    @staticmethod
    def split_feature_index_string(index_string: str):
        # splits string with form "45..68" into (45, 68)
        indices = index_string.strip().split("..")
        indices = [int(num) for num in indices]
        return tuple(indices)

    @staticmethod
    def split_closing_pair_index_string(index_string: str):
        # splits string with form "(45,68)" into (45, 68)
        indices = index_string.strip("(").strip(")").split(",")
        indices = [int(num) for num in indices]
        return tuple(indices)

    @staticmethod
    def split_closing_pair_base_string(base_string: str):
        # splits string with form "45:68" into (45, 68)
        bases = base_string.strip().split(":")
        return tuple(bases)

    @staticmethod
    def get_feature_sequence(feature_data: str) -> str:
        sequence = ''
        for char in feature_data:
            if char.isalpha():
                sequence += char

        return sequence

    @staticmethod
    def get_pseudoknot_label(feature_data: str):
        # get pk info
        if feature_data != '':
            return feature_data[3]
        return None

    @staticmethod
    def get_loop_parent_label(feature_data: str) -> str:
        parentLabel = ''
        for char in feature_data:
            if char == '.':
                break
            else:
                parentLabel += char

        return parentLabel

    @staticmethod
    def _create_stem_from_file_line(file_line: str) -> Stem:
        stem_data = file_line.strip().split(" ")

        stemLabel = stem_data[0]  # get stem identifier

        part5p_start, part5p_stop = StructureCreator.split_feature_index_string(
            stem_data[1])  # get stem 5p indices

        part5p_seq = StructureCreator.get_feature_sequence(
            stem_data[2])  # stem 5p sequence

        part3p_start, part3p_stop = StructureCreator.split_feature_index_string(
            stem_data[3])  # get stem 3p indices

        part3p_seq = StructureCreator.get_feature_sequence(
            stem_data[4])  # stem 3p sequence

        return Stem(stemLabel, part5p_seq, part3p_seq,
                    (part5p_start, part5p_stop), (part3p_start, part3p_stop))

    @staticmethod
    def _create_hairpin_from_file_line(file_line: str) -> Hairpin:
        hairpinData = file_line.strip().split(" ")

        # get hairpin label
        hairpinLabel = hairpinData[0]

        # get index of start of hairpin
        hairpin_indices = StructureCreator.split_feature_index_string(
            hairpinData[1])

        # get sequence of hairpin
        hairpin_seq = StructureCreator.get_feature_sequence(hairpinData[2])

        # get indices for closing base
        close_5_prime_index, close_3_prime_index = \
            StructureCreator.split_closing_pair_index_string(hairpinData[3])

        close_5_prime_base, close_3_prime_base = \
            StructureCreator.split_closing_pair_base_string(hairpinData[4])

        return Hairpin(hairpinLabel, hairpin_seq, hairpin_indices,
                       (close_5_prime_base, close_3_prime_base),
                       (close_5_prime_index, close_3_prime_index))

    @staticmethod
    def _create_bulge_from_file_line(file_line: str):
        bulgeData = file_line.strip().split(" ")

        # get bulge label
        bulgeLabel = bulgeData[0]

        # get start index of bulge
        bulge_start, bulge_stop = StructureCreator.split_feature_index_string(
            bulgeData[1])

        # get bulge sequence
        bulge_seq = StructureCreator.get_feature_sequence(bulgeData[2])

        # get index of 5' base in preceding pair
        precedingPair5pIndex, precedingPair3pIndex = \
            StructureCreator.split_closing_pair_index_string(bulgeData[3])

        # get 5' base of preceding pair
        precedingPair5pBase, precedingPair3pBase = \
            StructureCreator.split_closing_pair_base_string(bulgeData[4])

        # get index of 5' base in trailing pair
        trailingPair5pIndex, trailingPair3pIndex = \
            StructureCreator.split_closing_pair_index_string(bulgeData[5])

        # get 5' and 3' bases of trailing pair
        trailingPair5pBase, trailingPair3pBase = \
            StructureCreator.split_closing_pair_base_string(bulgeData[6])

        # need to make sure trailing pair is ordered correctly
        if(abs(bulge_stop - trailingPair5pIndex) >
                abs(bulge_stop - trailingPair3pIndex)):
            trailingBasePair = (trailingPair3pBase, trailingPair5pBase)
            trailingBasePairIndex = (trailingPair3pIndex, trailingPair5pBase)
        else:
            trailingBasePair = (trailingPair5pBase, trailingPair3pBase)
            trailingBasePairIndex = (trailingPair5pIndex, trailingPair3pIndex)

        return Bulge(bulgeLabel, bulge_seq, (bulge_start, bulge_stop),
                     (precedingPair5pBase, precedingPair3pBase),
                     (precedingPair5pIndex, precedingPair3pIndex),
                     trailingBasePair, trailingBasePairIndex)

    @staticmethod
    def _create_internal_loop_from_file_line(loop_1_file_line: str,
                                             loop_2_file_line: str):
        loop1Data = loop_1_file_line.strip().split(" ")
        loop2Data = loop_2_file_line.strip().split(" ")

        # get internal loop parent label
        parent_label = StructureCreator.get_loop_parent_label(loop1Data[0])

        # get loop 1 start and stop indices
        loop_1_indices = StructureCreator.split_feature_index_string(
            loop1Data[1])

        # get loop 1 start and stop indices
        loop_2_indices = StructureCreator.split_feature_index_string(
            loop2Data[1])

        # get loop sequences
        loop_1_sequence = StructureCreator.get_feature_sequence(loop1Data[2])
        loop_2_sequence = StructureCreator.get_feature_sequence(loop2Data[2])

        # get indices for loop closing pairs
        loop_1_closing_pair_indices = \
            StructureCreator.split_closing_pair_index_string(loop1Data[3])

        loop_2_closing_pair_indices = \
            StructureCreator.split_closing_pair_index_string(loop2Data[3])

        loop_1_closing_pair = StructureCreator.split_closing_pair_base_string(
            loop1Data[4])
        loop_2_closing_pair = StructureCreator.split_closing_pair_base_string(
            loop2Data[4])

        closing_pairs = ((loop_1_closing_pair[0], loop_1_closing_pair[1]),
                         (loop_2_closing_pair[1], loop_2_closing_pair[0]))

        return InternalLoop(parent_label, "1", "2",
                            loop_1_sequence, loop_2_sequence, loop_1_indices,
                            loop_2_indices, closing_pairs,
                            loop_1_closing_pair_indices,
                            loop_2_closing_pair_indices)

    @staticmethod
    def _create_external_loop_from_file_line(file_line: str):
        externalLoopData = file_line.strip().split(" ")

        externalLoopLabel = externalLoopData[0]

        external_loop_indices = StructureCreator.split_feature_index_string(
            externalLoopData[1])

        external_loop_sequence = StructureCreator.get_feature_sequence(
            externalLoopData[2])

        closing_pair_5p_indices = StructureCreator.split_closing_pair_index_string(
            externalLoopData[3])

        closing_pair_5p = StructureCreator.split_closing_pair_base_string(
            externalLoopData[4])

        closing_pair_3p_indices = StructureCreator.split_closing_pair_index_string(
            externalLoopData[5])

        closing_pair_3p = StructureCreator.split_closing_pair_base_string(
            externalLoopData[6])

        return ExternalLoop(externalLoopLabel, external_loop_sequence,
                            external_loop_indices,
                            closing_pair_5p, closing_pair_5p_indices,
                            closing_pair_3p, closing_pair_3p_indices)

    @staticmethod
    def _create_end_from_file_line(file_line: str):
        end_data = file_line.strip().split(" ")

        # get label for the feature
        end_label = end_data[0]

        end_indices = StructureCreator.split_feature_index_string(end_data[1])

        end_sequence = StructureCreator.get_feature_sequence(end_data[2])

        return End(end_label, end_sequence, end_indices)

    @staticmethod
    def _create_NCBP_from_file_line(file_line: str):
        ncbpData = file_line.strip().split(" ")

        # get label for ncbp
        ncbpLabel = ncbpData[0]

        # get index of 5' base
        base1Span = int(ncbpData[1])

        # get 5' base
        base1 = ncbpData[2]

        # get index 3' base
        base2Span = int(ncbpData[3])

        # get 3' base
        base2 = ncbpData[4]

        # get location of NCBP in other secondary structure
        if ncbpData[5] == '':
            loc = None
        else:
            loc = ncbpData[5]

        return NCBP(ncbpLabel, (base1, base2), (base1Span, base2Span), loc)

    @staticmethod
    def _create_multiloop_from_subcomponent_list(parent_label: str,
                                                 subcomponents: list):

        subunitLabels = []  # list to store all subunit labels
        sequences = {}  # dictionary to store all sequences
        spans = {}  # dictionary to store all sequence spans
        closingPairs = {}  # dictionary to store closing pairs
        closingPairsSpan = {}  # dictionary to store closing pairs spans

        for loop in subcomponents:
            loop = loop.strip().split(" ")

            # get inner loop subunit label and append to subunits list
            subunitLabel = loop[0][-1]
            subunitLabels.append(subunitLabel)

            loop_indices = StructureCreator.split_feature_index_string(loop[1])
            spans[subunitLabel] = loop_indices

            loop_sequence = StructureCreator.get_feature_sequence(loop[2])
            sequences[subunitLabel] = loop_sequence

            closing_pair_5p_indices = StructureCreator.split_closing_pair_index_string(
                loop[3])

            closing_pair_5p = StructureCreator.split_closing_pair_base_string(
                loop[4])

            closing_pair_3p_indices = StructureCreator.split_closing_pair_index_string(
                loop[5])

            closing_pair_3p = StructureCreator.split_closing_pair_base_string(
                loop[6])

            # add closing pairs and closing pair spans to dictionaries
            closingPairs[subunitLabel] = (closing_pair_5p, closing_pair_3p)
            closingPairsSpan[subunitLabel] = (
                closing_pair_5p_indices, closing_pair_3p_indices)

        return Multiloop(parent_label, subunitLabels, sequences, spans,
                         closingPairs, closingPairsSpan)

    @staticmethod
    def is_stem(feature_string: str) -> bool:
        return feature_string[0] == 'S' and feature_string[1].isdigit()

    @classmethod
    def create(cls, filename=None):
        if filename is None:
            print("No file provided")
            return

        # check that file is valid structure type
        if filename[-3::] != '.st':
            raise StructureFileNotProvided

        # try to open the provided file + error handling
        try:
            f = open(filename, 'r')
        except OSError:  # error finding or opening file
            print("Error finding or Accessing file.")
            return
        except Exception:  # something weird happened
            print(
                'Something unexpected ocurred when accessing the file')
            return

        # create new structure object
        new_structure = Structure(filename=filename)

        # Variables to validate all features have been read
        sequenceRead = False
        dotBracketRead = False
        structureArrayRead = False
        varnaRead = False

        # iterate through all of the lines in the file
        for line in f:

            if line[0] == '#':
                # get name of RNA molecule
                if line[0:6] == '#Name:':
                    new_structure._name = line[6:].strip()

                # get length of the RNA sequence
                elif line[0:8] == '#Length:':
                    new_structure._length = int(line[8:].strip().strip(','))
                    new_structure._componentArray = np.empty(
                        new_structure._length, dtype=object)

                # get page number for molecule
                elif line[0:12] == '#PageNumber:':
                    new_structure._pageNum = int(line[12:].strip())

                else:
                    continue

            # get actual RNA sequence
            elif (sequenceRead is False):
                new_structure._sequence = line.strip()
                sequenceRead = True

            # get Dot Bracket Notation for the molecule
            elif (dotBracketRead is False):
                new_structure._DBN = line.strip()
                dotBracketRead = True

            # get Annotated symbol form of the molecule
            elif (structureArrayRead is False):
                new_structure._structureArray = line.strip()
                structureArrayRead = True

            # varna notation for the molecule
            elif (varnaRead is False):
                new_structure._varna = line.strip()
                varnaRead = True

                if (sequenceRead and dotBracketRead and structureArrayRead
                        and varnaRead):
                    # when all identifying data has been parsed, read rest
                    # of the file and parse the StructureComponents
                    features = f.read()
                    break
                else:
                    print(f'File: {filename} is not proper .st format')
                    new_structure.resetStructure()

        # split rest of file contents into a list of strings
        features = features.split('\n')[:-1]
        i = 0
        while i < (len(features)):  # iterate through the individual string

            # stems
            if cls.is_stem(features[i]):
                new_stem = cls._create_stem_from_file_line(features[i])
                new_structure.addStem(new_stem.label(), new_stem)
                new_structure._addStemToComponentArray(new_stem)

            # Hairpins
            elif features[i][0] == 'H':
                new_hairpin = cls._create_hairpin_from_file_line(features[i])
                new_structure.addHairpin(new_hairpin.label(), new_hairpin)
                new_structure._addHairpinToComponentArray(new_hairpin)

            # Bulges
            elif features[i][0] == 'B':
                new_bulge = cls._create_bulge_from_file_line(features[i])
                new_structure.addBulge(new_bulge.label(), new_bulge)
                new_structure._addBulgeToComponentArray(new_bulge)

            # Inner Loops
            elif features[i][0] == 'I' and re.search('I\d{1,3}.1', features[i]):
                # pass both inner loop components
                new_internal_loop = cls._create_internal_loop_from_file_line(
                    features[i], features[i+1])
                new_structure.addInternalLoop(
                    new_internal_loop.label(), new_internal_loop)
                new_structure._addInternalLoopToComponentArray(
                    new_internal_loop)

            # MultiLoops
            elif features[i][0] == 'M':
                parentLabel = cls.get_loop_parent_label(
                    features[i])  # get parent label of the multiloop
                subcomponents = []  # array store multiloop subcomponents

                # linear probe for other multiloop subcomponents
                while cls.get_loop_parent_label(features[i]) == parentLabel:
                    subcomponents.append(features[i])
                    i += 1

                # parse the entire multiloop subcomponent list
                new_multiloop = cls._create_multiloop_from_subcomponent_list(
                    parentLabel, subcomponents)
                new_structure.addMultiLoop(
                    new_multiloop.label(), new_multiloop)
                new_structure._addMultiLoopToComponentArray(new_multiloop)
                continue

            # external loops
            elif features[i][0] == 'X':
                new_external_loop = cls._create_external_loop_from_file_line(
                    features[i])
                new_structure.addExternalLoop(
                    new_external_loop.label(), new_external_loop)
                new_structure._addExternalLoopToComponentArray(
                    new_external_loop)

            # NCBP
            elif features[i][0:4] == 'NCBP':
                new_NCBP = cls._create_NCBP_from_file_line(features[i])
                new_structure.addNCBP(new_NCBP.label(), new_NCBP)

            # Ends
            elif features[i][0] == 'E':
                new_end = cls._create_end_from_file_line(features[i])
                new_structure.addEnd(new_end.label(), new_end)
                new_structure._addEndToComponentArray(new_end)

            i += 1  # increment counter

        f.close()  # close the file

        new_structure._addStemBulgeNeighborBooleans()
        return new_structure
