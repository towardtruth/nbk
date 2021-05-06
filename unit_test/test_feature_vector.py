import unittest
import settings
from Data import FeatureVector
from models import Configuration, FeatureInfo
import sqls
from utils.DBManager import Pgsql


class TestFeatureVector(unittest.TestCase):
    # set target features
    #target_features = '1'
    #target_features = '23'
    #target_features = '1,2,3,4'
    #target_features = '4,7,9,10'
    #target_features = '1,2,3,4,5,6,7,8,9,10'
    #target_features = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22'
    #target_features = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20'
    target_features = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23'
    # set assigned genes limit
    #assigned_genes_limit = [int(10) + x for x in range(1, 24)]
    #assigned_genes_limit = [int(50) for x in range(1, 24)]
    assigned_genes_limit = [int(10) for x in range(1, 24)]
    exp_setting = Configuration()

    def test_feature_vector_random(self):
        self.exp_setting.set_target_features(target_features=self.target_features)

        self.exp_setting.set_assigned_genes_limit(self.assigned_genes_limit)

        # set debug_mode
        debug_mode = [1, 0]
        self.exp_setting.set_debug_mode(debug_mode)

        # set test_mode
        # exp_setting.set_test_mode(test_mode=True)
        # exp_setting.set_test_mode(test_mode=settings.TEST_MODE_KMER_FREQ)
        self.exp_setting.set_test_mode(test_mode=settings.TEST_MODE_TRUE)

        # set NEG class mode
        self.exp_setting.set_neg_class_mode(settings.NEG_CLASS_MODE_RND_S)

        # set class size
        self.exp_setting.set_class_size(class_size=2)

        # set fold size for cross-validation
        self.exp_setting.set_fold_size(fold_size=10)

        # set kmer size
        self.exp_setting.set_kmer_size(kmer_size=3)

        # set seq_type
        self.exp_setting.set_seq_type(seq_type='p')

        # set fsid
        fs_info = FeatureInfo(fsid=0,
                              fs_name='SM_RND',
                              gp_type='g',
                              class_size=2)
        self.exp_setting.set_fs_info(fs_info=fs_info)

        # Feature Vector
        fv = FeatureVector(self.exp_setting)

        # fv.test_features_dataset()
        fv.cross_validation()

        # write prediction results into a file (.csv)
        #fv.write_prediction_results()

        # Build feature vectors with prediction results
        #fv.build_feature_vector()
        print('Build Feature Vector')
        #gene_dataset = self.exp_setting.get_gene_dataset()
        all_gnids_list = self.exp_setting.get_gene_dataset_all_gnids_list()
        #for gnid in gene_dataset.get_all_gnids():
        for gnid in all_gnids_list:
            #for gnid in self.wd_all_gnids:
            gnid_vector = list()
            for feature in fv.features:
                #if feature.prediction_results is not None:
                if feature.prediction_results:
                    predicted_results = feature.prediction_results.get(gnid, None)
                    if predicted_results is None:
                        gnid_vector.append('?')
                    else:
                        #gnid_vector.append(feature.prediction_results[gnid].get_predicted_class())
                        gnid_vector.append(predicted_results.get_predicted_class())
            #print(gnid_vector)
            fv.feature_vector[gnid] = gnid_vector

        # TEST write feature vector
        for feature in fv.features:
            print('Feature Name: {}, tissue#: {}'.format(feature.name, feature.corresp_tissue))

            # gnids
            gnids = self.exp_setting.get_gene_dataset_gnids_list(feature_id=feature.corresp_tissue)

            #for gnid, vector in fv.feature_vector.items():
            for gnid in gnids:
                vector = fv.feature_vector.get(gnid, None)
                if vector is not None:
                    line = ",".join(str(value) for value in vector)
                    predicted_results = feature.prediction_results.get(gnid, None)
                    if predicted_results is None:
                        data_label = '?'
                    else:
                        data_label = predicted_results.get_assigned_class()
                    #f.write("%s,%s,%s\n" % (gnid, line, feature.prediction_results[gnid].get_assigned_class()))
                    print("%s,%s,%s\n" % (gnid, line, data_label))

                    #print('gnid:', gnid)
                    #print('line:', line)
                    #print('data_label:', data_label)

        #fv.write_feature_vector()
        #fv.create_arff()
        #fv.write_prediction_summary()

    def test_feature_vector_by_fsid(self):
        self.exp_setting.set_target_features(target_features=self.target_features)

        self.exp_setting.set_assigned_genes_limit(self.assigned_genes_limit)

        # set debug_mode
        debug_mode = [1, 0]
        self.exp_setting.set_debug_mode(debug_mode)

        # set test_mode
        # exp_setting.set_test_mode(test_mode=True)
        # exp_setting.set_test_mode(test_mode=settings.TEST_MODE_KMER_FREQ)
        self.exp_setting.set_test_mode(test_mode=settings.TEST_MODE_TRUE)

        # set NEG class mode
        self.exp_setting.set_neg_class_mode(settings.NEG_CLASS_MODE_RND_S)

        # set class size
        self.exp_setting.set_class_size(class_size=2)

        # set fold size for cross-validation
        self.exp_setting.set_fold_size(fold_size=10)

        # set kmer size
        self.exp_setting.set_kmer_size(kmer_size=3)

        # set seq_type
        self.exp_setting.set_seq_type(seq_type='p')

        # set fsid
        fsid = 44
        res_fs_info = Pgsql.Common.select_data(sqls.get_feature_set, (fsid))
        fs_info = FeatureInfo(fsid=fsid,
                              fs_name=res_fs_info[0][0].strip(),
                              gp_type=res_fs_info[0][1].strip(),
                              class_size=int(res_fs_info[0][2]))
        self.exp_setting.set_fs_info(fs_info=fs_info)

        # Feature Vector
        fv = FeatureVector(self.exp_setting)

        # fv.test_features_dataset()
        fv.cross_validation()

        # write prediction results into a file (.csv)
        #fv.write_prediction_results()

        # Build feature vectors with prediction results
        #fv.build_feature_vector()
        print('Build Feature Vector')
        #gene_dataset = self.exp_setting.get_gene_dataset()
        all_gnids_list = self.exp_setting.get_gene_dataset_all_gnids_list()
        #for gnid in gene_dataset.get_all_gnids():
        for gnid in all_gnids_list:
        #for gnid in self.wd_all_gnids:
            gnid_vector = list()
            for feature in fv.features:
                #if feature.prediction_results is not None:
                if feature.prediction_results:
                    predicted_results = feature.prediction_results.get(gnid, None)
                    if predicted_results is None:
                        gnid_vector.append('?')
                    else:
                        #gnid_vector.append(feature.prediction_results[gnid].get_predicted_class())
                        gnid_vector.append(predicted_results.get_predicted_class())
            #print(gnid_vector)
            fv.feature_vector[gnid] = gnid_vector

        # TEST write feature vector
        for feature in fv.features:
            print('Feature Name: {}, tissue#: {}'.format(feature.name, feature.corresp_tissue))

            # gnids
            gnids = self.exp_setting.get_gene_dataset_gnids_list(feature_id=feature.corresp_tissue)

            #for gnid, vector in fv.feature_vector.items():
            for gnid in gnids:
                vector = fv.feature_vector.get(gnid, None)
                if vector is not None:
                    line = ",".join(str(value) for value in vector)
                    predicted_results = feature.prediction_results.get(gnid, None)
                    if predicted_results is None:
                        data_label = '?'
                    else:
                        data_label = predicted_results.get_assigned_class()
                    #f.write("%s,%s,%s\n" % (gnid, line, feature.prediction_results[gnid].get_assigned_class()))
                    print("%s,%s,%s\n" % (gnid, line, data_label))

                    #print('gnid:', gnid)
                    #print('line:', line)
                    #print('data_label:', data_label)

        #fv.write_feature_vector()
        #fv.create_arff()
        #fv.write_prediction_summary()


if __name__ == '__main__':
    unittest.main()
