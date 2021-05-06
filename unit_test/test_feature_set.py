import unittest
from GeneGroups import FeatureSet


class TestFeatureSet(unittest.TestCase):
    def test_init_instance(self):
        fs = FeatureSet()
        print(fs.name)
        print(fs.fsid)
        print(fs.label_data_type)
        print(fs.description)
        print(fs.class_size)

    def test_get_feature_set(self):
        fs = FeatureSet.get(62)
        print(fs.name)
        print(fs.fsid)
        print(fs.label_data_type)
        print(fs.description)
        print(fs.class_size)


if __name__ == '__main__':
    unittest.main()
