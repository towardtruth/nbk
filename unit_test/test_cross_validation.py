import unittest
import settings
from Validation import CrossValidation
from GeneGroups import Genes
from models import Configuration


class TestCrossValidation(unittest.TestCase):
    gnids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    assigned_gnids = [1, 2, 3, 4, 5]
    corresp_tissue = 1
    genes = Genes(gnids)
    class_size = 2
    fold_size = 2
    kmer_size = 3
    target_features = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23'
    exp_setting = Configuration()
    exp_setting.set_target_features(target_features=target_features)
    #exp_setting.set_neg_class_mode(settings.NEG_CLASS_MODE_RND_S)
    exp_setting.set_neg_class_mode(settings.NEG_CLASS_MODE_NOT_P)

    def test_cross_validation(self):
        all_prediction_results = list()
        for corresp_tissue in range(1, 24):
            cv = CrossValidation(genes=self.genes,
                                 all_gnids=self.gnids,
                                 class_size=self.class_size,
                                 fold_size=self.fold_size,
                                 kmer_size=self.kmer_size,
                                 exp_setting=self.exp_setting)
            cv.build_datasets(assigned_genes=self.assigned_gnids,
                              neg_class_mode=self.exp_setting.get_neg_class_mode(),
                              corresp_tissue=corresp_tissue)
            prediction_results = cv.validation()
            all_prediction_results.append(prediction_results)

        feature_vector = dict()
        for gnid in self.gnids:
            gnid_vector = list()
            for prediction_results in all_prediction_results:
                if prediction_results:
                    gnid_pr = prediction_results.get(gnid)
                    if gnid_pr is None:
                        gnid_vector.append('?')
                    else:
                        gnid_vector.append(gnid_pr.get_predicted_class())
            feature_vector[gnid] = gnid_vector



if __name__ == '__main__':
    unittest.main()
