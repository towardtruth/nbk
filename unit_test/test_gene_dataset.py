import unittest
from nbk.models import GeneList, GeneDataSet, Configuration
import settings

# test list
list_a = [11, 3, 4, 5, 7, 8, 9, 1]
list_ax = [12, 3, 4, 5, 7, 8, 9, 1]
list_an = [12, 13, 14, 15, 17, 18, 19, 16]
list_b = [12, 14, 5, 6, 8, 10, 11, 2, 4]
list_bx = [2, 14, 5, 6, 8, 10, 11, 4]
list_bn = [112, 114, 15, 16, 18, 110, 111, 21, 41]
list_ab = [10, 12, 3, 14, 5, 6, 7, 8, 9, 1, 11, 2, 4]
list_c = [11,22,33,44,55]
list_cn = [77,23,38,45,79]
list_d = [66,77,88,99,100]
list_dn = [98,73,89,39,108]

class TestGeneDataset(unittest.TestCase):

    def test_gene_list(self):
        gene_list_a = GeneList(list_a)
        gene_list_b = GeneList(list_b)
        gene_list_ab = GeneList(list_ab)
        gene_list_merged = gene_list_a.merge_to(gene_list_b)
        print(gene_list_a.get_gene_list())
        print(gene_list_b.get_gene_list())
        print(gene_list_ab.get_gene_list())
        print(gene_list_merged.get_gene_list())
        self.assertEqual(list_a, list_a)
        self.assertNotEqual(list_a, list_ax)
        self.assertEqual(gene_list_ab.get_gene_list(), gene_list_merged.get_gene_list())

    def test_gene_dataset_same(self):
        # GeneDataSet TEST
        gene_dataset1 = GeneDataSet(positive=list_a, negative=list_b)
        gene_dataset2 = GeneDataSet(positive=list_a, negative=list_b)
        gene_list_ds1_pos = gene_dataset1.get_positive_class()
        gene_list_ds1_neg = gene_dataset1.get_negative_class()
        gene_list_ds1_all = gene_dataset1.get_all_genes()
        gene_list_ds2_pos = gene_dataset2.get_positive_class()
        gene_list_ds2_neg = gene_dataset2.get_negative_class()
        gene_list_ds2_all = gene_dataset2.get_all_genes()
        print('gene_list (ds1.pos): {}'.format(gene_list_ds1_pos.get_gene_list()))
        print('gene_list (ds1.neg): {}'.format(gene_list_ds1_neg.get_gene_list()))
        print('gene_list (ds2.pos): {}'.format(gene_list_ds2_pos.get_gene_list()))
        print('gene_list (ds2.neg): {}'.format(gene_list_ds2_neg.get_gene_list()))
        print('gene_list (ds1.all): {}'.format(gene_list_ds1_all.get_gene_list()))
        print('gene_list (ds2.all): {}'.format(gene_list_ds2_all.get_gene_list()))
        self.assertEqual(gene_list_ds1_pos.get_gene_list(), gene_list_ds2_pos.get_gene_list())
        self.assertEqual(gene_list_ds1_neg.get_gene_list(), gene_list_ds2_neg.get_gene_list())
        self.assertEqual(gene_list_ds1_all.get_gene_list(), gene_list_ds2_all.get_gene_list())

    def test_gene_dataset_not_same(self):
        # GeneDataSet TEST
        gene_dataset1 = GeneDataSet(positive=list_a, negative=list_b)
        gene_dataset2 = GeneDataSet(positive=list_ax, negative=list_bx)
        gene_list_ds1_pos = gene_dataset1.get_positive_class()
        gene_list_ds1_neg = gene_dataset1.get_negative_class()
        gene_list_ds1_all = gene_dataset1.get_all_genes()
        gene_list_ds2_pos = gene_dataset2.get_positive_class()
        gene_list_ds2_neg = gene_dataset2.get_negative_class()
        gene_list_ds2_all = gene_dataset2.get_all_genes()
        print('gene_list (ds1.pos): {}'.format(gene_list_ds1_pos.get_gene_list()))
        print('gene_list (ds1.neg): {}'.format(gene_list_ds1_neg.get_gene_list()))
        print('gene_list (ds2.pos): {}'.format(gene_list_ds2_pos.get_gene_list()))
        print('gene_list (ds2.neg): {}'.format(gene_list_ds2_neg.get_gene_list()))
        print('gene_list (ds1.all): {}'.format(gene_list_ds1_all.get_gene_list()))
        print('gene_list (ds2.all): {}'.format(gene_list_ds2_all.get_gene_list()))
        self.assertNotEqual(gene_list_ds1_pos.get_gene_list(), gene_list_ds2_pos.get_gene_list())
        self.assertNotEqual(gene_list_ds1_neg.get_gene_list(), gene_list_ds2_neg.get_gene_list())
        self.assertEqual(gene_list_ds1_all.get_gene_list(), gene_list_ds2_all.get_gene_list())

    def test_configuration_gene_dataset(self):
        exp_setting = Configuration()
        exp_setting.set_gene_dataset(feature_id=1, positive_class=list_a)
        exp_setting.set_gene_dataset(feature_id=1, negative_class=list_an)
        exp_setting.set_gene_dataset(feature_id=2, positive_class=list_b)
        exp_setting.set_gene_dataset(feature_id=2, negative_class=list_bn)
        exp_setting.set_gene_dataset(feature_id=3, positive_class=list_c)
        exp_setting.set_gene_dataset(feature_id=3, negative_class=list_cn)
        exp_setting.set_gene_dataset(feature_id=4, positive_class=list_d)
        exp_setting.set_gene_dataset(feature_id=4, negative_class=list_dn)

        # test - get a single dataset per feature_id
        gene_dataset_1 = exp_setting.get_gene_dataset(feature_id=1)
        gene_dataset_2 = exp_setting.get_gene_dataset(feature_id=2)
        gene_dataset_3 = exp_setting.get_gene_dataset(feature_id=3)
        gene_dataset_4 = exp_setting.get_gene_dataset(feature_id=4)
        print(gene_dataset_1.get_all_gnids())
        print(gene_dataset_2.get_all_gnids())
        print(gene_dataset_3.get_all_gnids())
        print(gene_dataset_4.get_all_gnids())

        # test - get all gnids in dataset list
        print(exp_setting.get_gene_dataset_all_gnids_list())

    def test_configuration_missing_gnids_in_promoter(self):
        exp_setting = Configuration()
        settings.conn_string = settings.conn_string_test
        missing_gnids = exp_setting.get_missing_gnids_in_promoter()
        print(missing_gnids)


if __name__ == '__main__':
    unittest.main()


