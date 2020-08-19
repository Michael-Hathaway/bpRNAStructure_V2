import unittest
import bpRNAStructure.StructureCreator as SC


class StructureParserUnittest(unittest.TestCase):
    def setUp(self):
        self.test_creator = SC.StructureCreator()

    def test_split_feature_index_string(self):
        test_string = "56..89"
        result = self.test_creator.split_feature_index_string(test_string)

        self.assertEqual(
            result, (56, 89), "Error in function: split_feature_index_string")

    def test_split_closing_pair_index_string(self):
        test_string = "(1089,1022)"
        result = self.test_creator.split_closing_pair_index_string(test_string)

        self.assertEqual(result, (1089, 1022),
                         "Error in function: split_closing_pair_index_string")

    def test_split_closing_pair_base_string(self):
        test_string = "A:G"
        result = self.test_creator.split_closing_pair_base_string(test_string)

        self.assertEqual(result, ("A", "G"),
                         "Error in function: split_closing_pair_base_string")


if __name__ == "__main__":
    unittest.main(verbosity=2)
