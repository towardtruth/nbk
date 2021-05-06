'''
    Unit Test - Build dataset with k-mer frequency features
    Author: Kyoung Tak Cho
    Created:
    Updated:

    Feature Sets:
      44  GE  N95
      45  GE  N90
      46  GE  N85
      47  GE  N80
      48  GE  N75
      49  GE  N70
'''
import unittest
from GeneGroups import Features, FeatureSet
from models import FreqDict
import sqls
from utils.DBManager import Pgsql
import pandas as pd
import random


class TestKmerFreqDataSet(unittest.TestCase):
    amino_acid = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    amino_acid_sorted = sorted(amino_acid)

    def _init_kmer_dict(self):
        kmer_dict = FreqDict(3)
        for i in self.amino_acid_sorted:
            for j in self.amino_acid_sorted:
                for k in self.amino_acid_sorted:
                    kmer = i+j+k
                    kmer_dict.add_kmer_freq(kmer, 0)
                    #print(kmer)
        # for test
        #print(len(kmer_dict))
        return kmer_dict

    def test_build_dataset_by_gnid(self):
        # get gnids per feature set number and tissue number
        seq_type = 'p'
        feature_set_id = 49
        tissue_id = 1
        feature_set = FeatureSet.get(feature_set_id)
        train_data = list()
        train_label = list()

        # positive class
        gnids_pos = Features.get_gnids_pos(feature_set_id, tissue_id)
        # negative class
        gnids_neg = Features.get_gnids_neg(feature_set_id=feature_set_id, seq_type=seq_type,
                                           label_type=feature_set.label_data_type,
                                           tissue_id=tissue_id,
                                           exclusive_gnids=gnids_pos)
        gnids = gnids_pos + gnids_neg
        print("pos: {}, neg: {}, total: {}".format(len(gnids_pos), len(gnids_neg), len(gnids)))

        #print(gnids, len(gnids))
        #for gnid in sorted(gnids):
        for gnid in gnids_pos:
            train_data.append(self._get_kmer_freq_by_gnid(gnid))
            #print('{} {}/{}'.format(gnid, len(train_data), len(gnids)))
            train_label.append(1)
        for gnid in gnids_neg:
            train_data.append(self._get_kmer_freq_by_gnid(gnid))
            #print('{} {}/{}'.format(gnid, len(train_data), len(gnids)))
            train_label.append(0)


        print(len(train_data))

        # shuffle data
        map_idx_pos = list(zip(train_data, train_label))
        random.shuffle(map_idx_pos)
        train_data, train_label = zip(*map_idx_pos)

        x_train = pd.DataFrame(train_data)
        # write x_train into CSV file
        x_train.to_csv('x_train.csv', index=False)

        y_train = pd.DataFrame(train_label)
        y_train.to_csv('y_train.csv', index=False)


    def _get_crsp_tissues(self):
        feature_set_id = 49
        tids = FeatureSet.tissues(feature_set_id)
        for tissue_id in tids:
            assigned_gnids = Features.get_gnids_pos(feature_set_id, tissue_id)
            print(assigned_gnids)

    def _get_kmer_freq_by_gnid(self, gnid=None):
        #if gnid is None:
        #    raise ValueError('gnid is empty.')

        kmer_dict = self._init_kmer_dict()
        #gnid = 1000
        seq_type = 'p'
        k = 3

        # get kmer frequency by gnid
        sql = sqls.get_kmer_freq_by_gnid
        pars = (gnid, seq_type, k)
        res = Pgsql.Common.select_data(sql=sql, pars=pars)
        gsid = None

        for row in res:
            kmer = str(row[0]).upper()
            freq = row[1]
            kmer_dict.add_kmer_freq(kmer, freq)
            if gsid is None:
                gsid = row[2]

        print("gnid: {}, gsid: {}".format(gnid, gsid))
        #kmer_dict.print(sort_type=1)
        freq_list_sorted = kmer_dict.to_list(key_value=1, sort_type=1)
        #print(freq_list_sorted)
        #return kmer_dict
        return freq_list_sorted


if __name__ == '__main__':
    unittest.main()
