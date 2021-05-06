"""
    UNIT TEST - k-mer distribution
    Author: Kyoung Tak Cho
"""

import unittest
from nbk.GeneGroups import Features, FeatureSet
from nbk import sqls
from nbk import settings
from nbk.utils.DBManager import Pgsql
from nbk.utils.FileManager import ArffManager
from nbk.utils.Common import ListTool
from nbk.Kmer import Kmer
#import numpy as np
from km_distance.KmerDistance import PositionVector, KmerPair, KmerDistance
from nbk.Data import GeneSeqContainer, DataSetContainer


#
# SQLs
#


class TestKmerDistribution(unittest.TestCase):
    def test_get_gnids(self):
        # -. get gnids by tissue from feature table
        assigned_gnids = Features.get_gnids_pos(62, 1)
        print(assigned_gnids)

    def test_get_crsp_tissues(self):
        tids = FeatureSet.tissues(62)
        for tissue_id in tids:
            assigned_gnids = Features.get_gnids_pos(62, tissue_id)
            print(assigned_gnids)

    def test_get_feature_set(self):
        fsids = [56, 57, 58, 59, 60, 61]
        for fsid in fsids:
            fs = FeatureSet.get(fsid)
            #print(fs.name)
            #print(fs.fsid)
            #print(fs.label_data_type)
            #print(fs.class_size)
            print(fs.exp_level())
            print(fs.exp_type())

    def test_build_dataset(self):
        use_max_seq_len = True
        seq_type = 'p'
        #fsids = [56, 57, 58, 59, 60, 61]
        fsids = [62, 63, 64, 65, 66, 67,
                 68, 69, 70, 71, 72, 73,
                 74, 75, 76, 77, 78, 79,
                 99, 100, 101, 102, 103,
                 104, 105, 106, 107, 108, 109, 110]

        dataset = list()
        for fsid in fsids:
            fs = FeatureSet.get(fsid)
            tids = FeatureSet.tissues(fsid)

            # exp_type
            exp_type = None
            if fs.exp_type().rstrip() == "g":
                exp_type = "RA"
            elif fs.exp_type().rstrip() == "p":
                exp_type = "PA"
            elif fs.exp_type().rstrip() == "b":
                exp_type = "CA"

            # exp_level
            exp_level = fs.exp_level().rstrip().replace('N', 'T')

            print(exp_type)
            print(exp_level)

            conn = Pgsql.Common.connect()
            conn_kmd = Pgsql.Common.connect(settings.conn_string_kmd)
            for tissue_id in tids:
                # get gnids for positive class
                assigned_gnids = sorted(Features.get_gnids_pos(fsid, tissue_id, cur=conn.cursor()))
                #print(assigned_gnids)

                # get gsids for positive class
                pos_class = sorted(Features.get_gsids_pos(assigned_gnids, seq_type, use_max_seq_len, cur=conn.cursor()))
                #print(pos_class)

                # get gnids for negative class
                random_gnids = sorted(Features.get_gnids_neg(feature_set_id=fsid, seq_type=seq_type,
                                                             label_type=fs.exp_type().rstrip(), tissue_id=tissue_id,
                                                             exclusive_gnids=assigned_gnids))
                #print(random_gnids)

                # get gsids for negative class
                neg_class = sorted(Features.get_gsids_pos(random_gnids, seq_type, use_max_seq_len, cur=conn.cursor()))
                #print(neg_class)

                # sql for insert dataset table for pos&neg genes
                #  - tissue, exp_type, exp_level, name, pos, neg, desc

                # dataset name
                name = "{:02d}-{}-{}".format(tissue_id, exp_type, exp_level)
                pos_class_str = ListTool.list2str(pos_class)
                neg_class_str = ListTool.list2str(neg_class)
                #dataset.append((tissue_id, exp_type, exp_level, name, pos_class_str, neg_class_str))
                par = (tissue_id, exp_type, exp_level, name, pos_class_str, neg_class_str)
                Pgsql.Common.insert_data(sql=sqls.build_dataset, par_list=par, conn=conn_kmd)
                print(name)

            conn.close()
            conn_kmd.close()

        #for ds in dataset:
        #    print(ds)
        #    Pgsql.Common.insert_data(sql=sqls.build_dataset, par_list=ds)
        print(len(dataset))

    def test_build_kmer_freq(self):
        # get gsids per dataset
        conn_kmd = Pgsql.Common.connect(settings.conn_string_kmd)

        #exp_types = "'RA'"
        exp_types = "'PA'"
        #exp_types = "'RA', 'PA'"
        exp_levels = "'T95', 'T90', 'T85', 'T80', 'T75', 'T70'"
        tissues = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23'
        #tissues = '1'#,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23'
        k = 4
        cnt_test = 0

        res = Pgsql.Common.select_data(sql=sqls.get_dataset,
                                       pars=(exp_types, exp_levels, tissues), cur=conn_kmd.cursor())
        for row in res:
            dsid = row[0]
            name = row[4]
            pos_class_str = row[5]
            neg_class_str = row[6]
            #pos_class = list(map(int, pos_class_str.split(',')))
            #neg_class = list(map(int, neg_class_str.split(',')))
            #print(len(pos_class))

            # positive class
            seqs_pos = Pgsql.Common.select_data(sql=sqls.get_seq_by_gsids,
                                                pars=(pos_class_str), cur=conn_kmd.cursor())
            kmer_freq = None
            for seq in seqs_pos:
                kmer_freq = Kmer.kmer_freq_acc(seq=seq[2], window_size=k, freq=kmer_freq)

            total_len = len(kmer_freq)
            cnt = 1
            par_list = list()
            for kmer, freq in sorted(kmer_freq.items(), key=lambda kv: (kv[1], kv[0]), reverse=True):
                percentile = (total_len - cnt + 1) / total_len * 100.0
                par_list.append((dsid, 1, percentile, freq, kmer, k))
                #print(cnt, kmer, freq, percentile)
                cnt += 1

            # negative class
            seqs_neg = Pgsql.Common.select_data(sql=sqls.get_seq_by_gsids,
                                                pars=(neg_class_str), cur=conn_kmd.cursor())
            kmer_freq = None
            for seq in seqs_neg:
                kmer_freq = Kmer.kmer_freq_acc(seq=seq[2], window_size=k, freq=kmer_freq)

            total_len = len(kmer_freq)
            cnt = 1
            #par_list = list()
            for kmer, freq in sorted(kmer_freq.items(), key=lambda kv: (kv[1], kv[0]), reverse=True):
                percentile = (total_len - cnt + 1) / total_len * 100.0
                par_list.append((dsid, 0, percentile, freq, kmer, k))
                #print(cnt, kmer, freq, percentile)
                cnt += 1

            # insert kmer freq into the DB
            for par in par_list:
                Pgsql.Common.insert_data(sql=sqls.build_kmer_freq, par_list=par, conn=conn_kmd)
                #cnt_test += 1
            print("DSID: {}, {} :: {} rows inserted.".format(dsid, name.rstrip(), len(par_list)))
            #print(cnt_test)
        #print(cnt_test)

        conn_kmd.close()

    def test_data_set_container(self):
        tissues = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23'
        exp_types = "'RA'"
        exp_levels = "'T95', 'T90', 'T85', 'T80', 'T75', 'T70'"
        datasets = DataSetContainer(exp_types, exp_levels, tissues)
        for ds in datasets:
            print(ds.dsid, ds.tissue, ds.name)


    def test_kmer_select(self):
        # TODO: biuld k-mer selection module
        # 1. Build accumulated kmer frequency table
        # 2. pick some k-mers at random
        # 3. compute k-mer distance matrix
        # 4. build dataset
        # 5. apply ML models such as BN and deep learning

        conn_kmd = Pgsql.Common.connect(settings.conn_string_kmd)

        # get target dsids
        #tissues = '1'#,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23'
        tissues = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23'
        #exp_types = "'RA'"
        exp_types = "'PA'"
        #exp_types = "'RA', 'PA'"
        #exp_levels = "'T95'"#, 'T90', 'T85', 'T80', 'T75', 'T70'"
        #exp_levels = "'T95', 'T90'"#, 'T85', 'T80', 'T75', 'T70'"
        #exp_levels = "'T90', 'T85', 'T80', 'T75', 'T70'"
        exp_levels = "'T95', 'T90', 'T85', 'T80', 'T75', 'T70'"
        #k = (3, 4, 5)
        #limit = (10, 15, 20)
        k = 3
        limit = 30

        datasets = DataSetContainer(exp_types, exp_levels, tissues)
        #res = Pgsql.Common.select_data(sql=sqls.get_dataset,
        #                               pars=(exp_types, exp_levels, tissues), cur=conn_kmd.cursor())

        # get significant k-mers for each dsid
        #for row in res:
        for ds in datasets:
            distance_vectors = list()

            #dsid = row[0]
            #name = row[4].rstrip()
            #dsid = ds.dsid
            #name = ds.name

            # positive_class
            kmers_pos = Pgsql.Common.select_data(sql=sqls.get_kmer_by_top_rank,
                                                 pars=(ds.dsid, k, '1', limit),
                                                 cur=conn_kmd.cursor(), to_list=True)
            kmers_pos = [kmer.rstrip() for kmer in kmers_pos]
            #print(kmers_pos)

            # negative_class
            kmers_neg = Pgsql.Common.select_data(sql=sqls.get_kmer_by_top_rank,
                                                 pars=(ds.dsid, k, '0', limit),
                                                 cur=conn_kmd.cursor(), to_list=True)
            kmers_neg = [kmer.rstrip() for kmer in kmers_neg]
            #print(kmers_neg)


            #
            # compute k-mer pair matrix
            #
            gene_seqs_pos = GeneSeqContainer(ds.pos_genes)
            gene_seqs_neg = GeneSeqContainer(ds.neg_genes)

            is_set_kp_label_pos = False
            is_set_kp_label_neg = False
            attributes = list()

            # Positive class
            for gene_seq in gene_seqs_pos:
                kd_pos = KmerDistance(gene_seq=gene_seq, km_set=kmers_pos)
                kd_neg = KmerDistance(gene_seq=gene_seq, km_set=kmers_neg)
                kd_pos._init_dict()
                kd_pos.distance_matrix()
                kd_neg._init_dict()
                kd_neg.distance_matrix()

                # show contents
                #kd.show_kp_set()
                #kd.show_kp_dict()
                #print(kd.km_distance)
                #print(len(kd.km_distance))
                #print("POS - ", kd_pos.distance_mean_vector())
                #print("NEG - ", kd_neg.distance_mean_vector())
                #print(kd_pos.distance_mean_vector() + kd_neg.distance_mean_vector())

                # build distance vector: [gnid, pos_distance_vector, neg_distance_vecotr, class]
                distance_vectors.append([gene_seq.gnid] + kd_pos.distance_mean_vector() + kd_neg.distance_mean_vector() + [1])

                # set attributes
                if not is_set_kp_label_pos:
                    #attributes += kd_pos.kmer_pair_labels()
                    attributes += ["{}_1".format(label) for label in kd_pos.kmer_pair_labels()]
                    #print(attributes)
                    is_set_kp_label_pos = True

            # Negative class
            for gene_seq in gene_seqs_neg:
                kd_pos = KmerDistance(gene_seq=gene_seq, km_set=kmers_pos)
                kd_neg = KmerDistance(gene_seq=gene_seq, km_set=kmers_neg)
                kd_pos._init_dict()
                kd_pos.distance_matrix()
                kd_neg._init_dict()
                kd_neg.distance_matrix()

                distance_vectors.append([gene_seq.gnid] + kd_pos.distance_mean_vector() + kd_neg.distance_mean_vector() + [0])

                # set attributes
                if not is_set_kp_label_neg:
                    attributes += ["{}_0".format(label) for label in kd_neg.kmer_pair_labels()]
                    #print(attributes)
                    is_set_kp_label_neg = True

            #for dv in distance_vectors:
            #    print(dv)
            #print(len(distance_vectors))

            #print(attributes)

            # save arff file
            d_name = "kmd_K{}_{}_S{:03d}".format(k, ds.get_name(type=1), limit)        # kmd_[dataset_name]_[k-size]_
            file_name = "../../kmd/results/{}.arff".format(d_name)
            ArffManager.write(file_name, d_name, attributes, distance_vectors)
            print("{} created.".format(file_name))

        conn_kmd.close()

    def test_get_kmer_distribution(self):
        kmer_dist = list()
        use_max_seq_len = True
        #seq_type = 'p'
        kmer_size = 3
        fsids = [62,63,64]

        for fsid in fsids:
            fs_info = FeatureSet.get(fsid)
            fs_name = fs_info.name
            seq_type = fs_info.label_data_type

            tids = FeatureSet.tissues(fsid)
            conn = Pgsql.Common.connect()
            for tissue_id in tids:
                assigned_gnids = Features.get_gnids_pos(fsid, tissue_id, cur=conn.cursor())

                # get gsids
                assigned_gsids = Features.get_gsids_pos(assigned_gnids, seq_type, use_max_seq_len, cur=conn.cursor())
                assigned_gsids = sorted(assigned_gsids)
                str_assigned_gsids = ','.join(str(gsid) for gsid in assigned_gsids)

                # get kmer distribution
                res = Pgsql.Common.select_data(sql=sqls.get_kmer_freq_sum,
                                                 pars=(kmer_size, str_assigned_gsids),
                                                 cur=conn.cursor())
                #print(res)
                kmer_dist.append([fsid, fs_name, tissue_id] + res)


            # write kmer_dist file (.csv)
            file_name = 'data_anl/{}_{}.csv'.format(fsid, fs_name)
            ListTool.list2csv(kmer_dist, file_name)
            print('{} saved.'.format(file_name))

            conn.close()

    def test_adding_kmers(self):
        pass

    def test_build_kmer_dictionary(self):
        pass



if __name__ == '__main__':
    unittest.main()

