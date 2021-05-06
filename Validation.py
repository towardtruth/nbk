'''
    Validation class
        1. K-fold Cross Validation
    Author: Kyoung Tak Cho
    Created: Tue Sep  4 04:03:20 CDT 2018
    Updated: Mon Jun 10 17:08:10 CDT 2019
'''

import settings
import sqls
import math
from models import ValidationDataset, FreqDict, GenePrediction, ConfusionMatrix
from GeneGroups import FreqDictSet, Genes
from utils.Common import ListTool
from utils.DBManager import Pgsql


#
# Validation - Kmer Frequency group for validation
#
class ValidationFreqDict(object):
    def __init__(self, datasets=None, k=None, genes=None):
        if datasets is None:
            error_mesg = 'datasets is empty.'
            raise ValueError(error_mesg)
        if k is None:
            error_mesg = 'k (kmer window size) is empty.'
            raise ValueError(error_mesg)

        self.datasets = datasets
        self.class_size = len(self.datasets)
        self.k = k
        if genes is None:
            self.genes = Genes()
        else:
            self.genes = genes
        self.freq_dict_set = list()     # FreqDictSet * class_size

        #
        # Initialization
        #

        # Build kmer frequency dictionary for each fold
        self._init_freq_dict_set()

    # Build kmer frequency dictionary for each fold
    def _init_freq_dict_set(self):
        for class_num, v_dataset in enumerate(self.datasets):
            print('class_num: ', class_num)
            self.freq_dict_set.append(list())
            self.freq_dict_set[class_num] = FreqDictSet()

            for fold_idx, fold in enumerate(v_dataset):
                print('fold_idx: ', fold_idx, 'len: ', len(fold))
                fd_k = FreqDict(k=self.k)
                fd_k_1 = FreqDict(k=self.k - 1)
                for gnid in fold:
                    gene_info = self.genes.get_gene(gnid=gnid)

                    # for k-mer
                    gene_fd_k = gene_info.pep_seq_max.get_kmer_freq(self.k)
                    for kmer, freq in gene_fd_k.kmer_freq.items():
                        fd_k.add_kmer_freq(kmer, freq)

                    # for (k-1)-mer
                    gene_fd_k_1 = gene_info.pep_seq_max.get_kmer_freq(self.k - 1)
                    for kmer, freq in gene_fd_k_1.kmer_freq.items():
                        fd_k_1.add_kmer_freq(kmer, freq)

                self.freq_dict_set[class_num].append_kmer_dict(k=self.k, freq_dict=fd_k)
                self.freq_dict_set[class_num].append_kmer_dict(k=self.k - 1, freq_dict=fd_k_1)

    def get_genes(self):
        return self.genes

    def get_fd_class(self, class_num=None, k=None):
        if class_num is None:
            error_mesg = 'class num is empty.'
            raise ValueError(error_mesg)
        if k is None:
            error_mesg = 'k is empty'
            raise ValueError(error_mesg)

        return self.freq_dict_set[class_num].get_fd_class(k=k)

    def get_fd_fold(self, class_num=None, k=None, fold_idx=None):
        if class_num is None:
                error_mesg = 'class num is empty.'
                raise ValueError(error_mesg)
        if k is None:
            error_mesg = 'k is empty'
            raise ValueError(error_mesg)
        if fold_idx is None:
            error_mesg = 'fold_idx is empty'
            raise ValueError(error_mesg)

        return self.freq_dict_set[class_num].get_fd_fold(k=k, fold_idx=fold_idx)

    def get_train_set(self, class_num=None, k=None, fold_idx=None):
        return self.freq_dict_set[class_num].get_fd_complement(k=k, fold_idx=fold_idx)

    def test_freq_dict_set(self):
        for class_num, v_dataset in enumerate(self.datasets):
            print('TEST - class_num: ', class_num)
            self.freq_dict_set.append(list())
            self.freq_dict_set[class_num] = FreqDictSet()

            for fold_idx, fold in enumerate(v_dataset):
                print('TEST - fold_idx: ', fold_idx, 'len: ', len(fold),)
                fd_k = FreqDict(k=self.k)
                fd_k_1 = FreqDict(k=self.k - 1)
                for gnid in fold:
                    gene_info = self.genes.get_gene(gnid=gnid)

                    # for k-mer
                    gene_fd_k = gene_info.pep_seq_max.get_kmer_freq(self.k)
                    for kmer, freq in gene_fd_k.kmer_freq.items():
                        fd_k.add_kmer_freq(kmer, freq)

                    # for (k-1)-mer
                    gene_fd_k_1 = gene_info.pep_seq_max.get_kmer_freq(self.k - 1)
                    for kmer, freq in gene_fd_k_1.kmer_freq.items():
                        fd_k_1.add_kmer_freq(kmer, freq)
                print('fd_k value sum: {}'.format(fd_k.total_freq()))
                print('fd_k_1 value sum: {}'.format(fd_k_1.total_freq()))

                self.freq_dict_set[class_num].append_kmer_dict(k=self.k, freq_dict=fd_k)
                self.freq_dict_set[class_num].append_kmer_dict(k=self.k - 1, freq_dict=fd_k_1)


#
# K-fold Cross Validation
#
class CrossValidation(object):
    def __init__(self, genes=None, all_gnids=None,
                 class_size=None, fold_size=None, kmer_size=None,
                 exp_setting=None):
        if exp_setting is not None and not exp_setting.get_test_mode():
            if genes is None:
                error_mesg = 'genes is empty.'
                raise ValueError(error_mesg)
            if class_size is None:
                error_mesg = 'class size is empty.'
                raise ValueError(error_mesg)
            if fold_size is None:
                error_mesg = 'fold size is empty.'
                raise ValueError(error_mesg)
            if kmer_size is None:
                error_mesg = 'kmer size is empty.'
                raise ValueError(error_mesg)

        # general attributes
        self.genes = genes              # Genes type
        self.all_gnids = all_gnids      # list() type
        self.datasets = list()          # ValidationDataset
        self.class_size = class_size
        self.fold_size = fold_size
        self.kmer_size = kmer_size
        self.prediction_results = dict()    # * number of gnids
        self.assigned_genes = list()
        self.exp_setting = exp_setting

    def build_datasets(self, assigned_genes=None, neg_class_mode=1, corresp_tissue=None):
        if assigned_genes is None:
            error_mesg = 'assigned is empty.'
            raise ValueError(error_mesg)
        else:   # build assigned gnids for promoter data
            if self.exp_setting.get_seq_type() == 'm1':
                # remove missing gnids for promoter data
                self.assigned_genes = self.remove_missing_gnids(assigned_genes)

                # set exclusive gnids
                exclusive_gnids = self.get_exclusive_gnids(assigned_genes)
            else:
                self.assigned_genes = assigned_genes.copy()
                exclusive_gnids = self.assigned_genes

        if corresp_tissue is None:
            raise ValueError('corresp_tissue is empty.')

        if self.fold_size is None:
            self.fold_size = 1
            print('fold size is empty. fold_size = 1 (default)')

        if self.fold_size <= 0:
            self.fold_size = 1
            print('fold size is less than zero. fold_size = 1 (default)')

        for class_num in range(0, self.class_size):
            if settings.NEED_NEW_GENE_DATA:
                self.datasets.append(class_num)
                self.datasets[class_num] = list()
                sub_dataset = list()
                gp_type = self.exp_setting.get_gp_type()
                if class_num == 0:
                    if neg_class_mode == settings.NEG_CLASS_MODE_NOT_P:
                        wd_all_gnids_per_tissue = GetData.wd_all_gnid_per_tissue(self.exp_setting, corresp_tissue)
                        sub_dataset = list(set(wd_all_gnids_per_tissue) - set(self.assigned_genes))

                    elif neg_class_mode == settings.NEG_CLASS_MODE_RND_S:   # random and same number with POS class
                        if gp_type not in ('g', 'p'):   # default is type is 'p' for gp_type == 'm1' (combined type)
                            gp_type = 'p'
                        pars = (gp_type, corresp_tissue, ListTool.list2str(exclusive_gnids, ','), len(self.assigned_genes))
                        random_assigned_gene = Pgsql.Common.select_data(sqls.get_gene_tissues_random, pars)
                        sub_dataset = ListTool.twoD2oneD(random_assigned_gene)

                    elif neg_class_mode == settings.NEG_CLASS_MODE_RND_M:
                        wd_all_gnids_per_tissue = GetData.wd_all_gnid_per_tissue(self.exp_setting, corresp_tissue)
                        sub_dataset = list(set(wd_all_gnids_per_tissue) - set(self.assigned_genes))

                    else:       # default
                        wd_all_gnids_per_tissue = GetData.wd_all_gnid_per_tissue(self.exp_setting, corresp_tissue)
                        sub_dataset = list(set(wd_all_gnids_per_tissue) - set(self.assigned_genes))

                    # set negative class genes
                    self.exp_setting.set_gene_dataset(feature_id=corresp_tissue, negative_class=sub_dataset)

                elif class_num == 1:
                    sub_dataset = self.assigned_genes

                    # set positive class genes
                    self.exp_setting.set_gene_dataset(feature_id=corresp_tissue, positive_class=sub_dataset)
                self.datasets[class_num] = ValidationDataset(gnids=sub_dataset, fold_size=self.fold_size)
            else:
                sub_dataset_neg = self.exp_setting.get_gene_dataset_neg()
                sub_dataset_pos = self.exp_setting.get_gene_dataset_pos()
                self.datasets[0] = ValidationDataset(gnids=sub_dataset_neg, fold_size=self.fold_size)
                self.datasets[1] = ValidationDataset(gnids=sub_dataset_pos, fold_size=self.fold_size)

    def remove_missing_gnids(self, all_gnids):
        # remove missing gnids for promoter
        missing_gnids = self.exp_setting.get_missing_gnids_in_promoter()
        #new_sub_dataset = list(set(all_gnids) - set(missing_gnids))
        new_sub_dataset = ListTool.sub(all_gnids, missing_gnids)
        print('all gnids:', len(all_gnids))
        print('common items:', len(ListTool.common_items(all_gnids, missing_gnids)))
        print('removed:', len(new_sub_dataset))
        return sorted(new_sub_dataset)

    def get_exclusive_gnids(self, all_gnids):
        missing_gnids = self.exp_setting.get_missing_gnids_in_promoter()
        exclusive_gnids = ListTool.add_rm_dup(all_gnids, missing_gnids)
        return sorted(exclusive_gnids)

    def test_datasets(self):
        if len(self.datasets) == 0:
            raise ValueError('dataset is empty.')

        for class_num, dataset in enumerate(self.datasets):
            for fold_idx, fold in enumerate(dataset):
                print('fold#: {}, genes: {}'.format(fold_idx, ','.join(str(gnid) for gnid in fold)))

    def validation(self):

        #
        # Phase 1: build Kmer Frequency dictionary for each fold
        #
        print('Step 1 - build frequency dictionaries')
        v_freq_dict = ValidationFreqDict(datasets=self.datasets, k=self.kmer_size, genes=self.genes)
        # update genes info
        self.genes = v_freq_dict.get_genes()

        # TEST
        if self.exp_setting.get_test_mode() == settings.TEST_MODE_KMER_FREQ:
            v_freq_dict.test_freq_dict_set()

        #
        # Phase 2: compute probabilities for each fold
        #
        print('Stepe 2 - prediction')

        # Gene prediction results dictionary
        prediction_results = dict()

        # build frequency dictionaries by class level
        fd_class_k = list()
        fd_class_k_1 = list()
        for class_num in range(0, self.class_size):
            fd_class_k.append(v_freq_dict.get_fd_class(class_num=class_num, k=self.kmer_size))
            fd_class_k_1.append(v_freq_dict.get_fd_class(class_num=class_num, k=self.kmer_size - 1))
            print('total: fd_class_k:', fd_class_k[class_num].total_freq_float(), len(fd_class_k[class_num]), len(fd_class_k))
            print('total: fd_class_k-1:', fd_class_k_1[class_num].total_freq_float(), len(fd_class_k_1[class_num]), len(fd_class_k_1))

        # set dataset size for each class and compute probability for each class
        p_class = [0., 0.]
        size_total = len(self.datasets[0]) + len(self.datasets[1])
        for class_num in range(0, self.class_size):
            p_class[class_num] = len(self.datasets[class_num]) / size_total

        for assigned_class, dataset in enumerate(self.datasets):

            for fold_idx, fold in enumerate(dataset):
                if self.exp_setting.get_test_mode() == settings.TEST_MODE_KMER_FREQ:
                    print('Total gnids in fold# {}: {}, gnids: {}'.format(fold_idx,
                                                                          len(fold),
                                                                          ','.join(str(gnid) for gnid in fold)))

                fd_train_set_k = list()     # * class_size
                fd_train_set_k_1 = list()   # * class_size
                fd_train_set_total_freq_k = list()
                fd_train_set_total_freq_k_1 = list()
                for class_num in range(0, self.class_size):
                    fd_train_set_k.append(fd_class_k[class_num] - v_freq_dict.get_fd_fold(class_num=class_num, k=self.kmer_size, fold_idx=fold_idx))
                    fd_train_set_k_1.append(fd_class_k_1[class_num] - v_freq_dict.get_fd_fold(class_num=class_num, k=self.kmer_size - 1, fold_idx=fold_idx))

                    fd_train_set_total_freq_k.append(fd_train_set_k[class_num].total_freq_float())
                    fd_train_set_total_freq_k_1.append(fd_train_set_k_1[class_num].total_freq_float())

                # initialize confusion matrix
                cm = ConfusionMatrix()

                # check for both classes
                if fd_train_set_total_freq_k[class_num] == 0 or fd_train_set_total_freq_k_1[class_num] == 0:
                    dataset.set_confusion_matrix(cm=cm, fold_idx=fold_idx)
                    prediction_results = dict()
                    return prediction_results

                for gnid in fold:

                    gene_info = self.genes.get_gene(gnid=gnid)
                    fd_gene_k = gene_info.pep_seq_max.get_kmer_freq(self.kmer_size)
                    fd_gene_k_1 = gene_info.pep_seq_max.get_kmer_freq(self.kmer_size - 1)
                    log_p = [0., 0.]

                    # 3-mer
                    for kmer in fd_gene_k.kmer_freq:
                        freq_k = fd_gene_k.kmer_freq[kmer]
                        if freq_k > 0:
                            for class_num in range(0, self.class_size):
                                freq_k_in_train = fd_train_set_k[class_num].get_kmer_freq_value(kmer=kmer)
                                freq_k_total = fd_train_set_total_freq_k[class_num]
                                prob = freq_k * math.log(freq_k_in_train / freq_k_total)
                                log_p[class_num] += prob

                                # TEST - KmerFreq/FreqDict
                                if self.exp_setting.get_test_mode() == settings.TEST_MODE_KMER_FREQ:
                                    print('class#: {}, kmer: {}, freq_k: {}, probability: {}, freq_k in train set: {}, freq_k total: {}'.format(
                                          class_num, kmer, freq_k, prob, freq_k_in_train, freq_k_total))
                        else:
                            print('kmer: {}, frequency: {}'.format(kmer, freq_k))
                            raise(ValueError('kmer frequency is less than 0'))

                    # 2-mer
                    for kmer in fd_gene_k_1.kmer_freq:
                        freq_k_1 = fd_gene_k_1.kmer_freq[kmer]
                        if freq_k_1 > 0:
                            for class_num in range(0, self.class_size):
                                freq_k_1_in_train = fd_train_set_k_1[class_num].get_kmer_freq_value(kmer=kmer)
                                freq_k_1_total = fd_train_set_total_freq_k_1[class_num]
                                prob = freq_k_1 * math.log(freq_k_1_in_train / freq_k_1_total)
                                log_p[class_num] -= prob

                                # TEST - KmerFreq/FreqDict
                                if self.exp_setting.get_test_mode() == settings.TEST_MODE_KMER_FREQ:
                                    print('class#: {}, kmer: {}, freq_k: {}, probability: {}, freq_k_1 in train set: {}, freq_k_1 total: {}'.format(
                                        class_num, kmer, freq_k_1, prob, freq_k_1_in_train, freq_k_1_total))

                    # log_p finalize - adding probability for each class
                    for class_num in range(0, self.class_size):
                        log_p[class_num] += math.log(p_class[class_num])

                    # TEST - KmerFreq/FreqDict
                    if self.exp_setting.get_test_mode() == settings.TEST_MODE_KMER_FREQ:
                        print('seq: {}'.format(gene_info.pep_seq_max.get_seq_str()))
                        print('log_p[0]: {}, log_p[1]: {}'. format(log_p[0], log_p[1]))

                    # set prediction results
                    gene_prediction = prediction_results.get(gnid, None)
                    if gene_prediction is None:
                        gene_prediction = GenePrediction(gnid=gnid, class_size=self.class_size, assigned_class=assigned_class)

                    if log_p[0] > log_p[1]:
                        predicted_class = 0

                        # set confusion_matrix
                        if assigned_class == 0:
                            cm.add_tn()
                        else:
                            cm.add_fn()
                    else:
                        predicted_class = 1

                        # set confusion_matrix
                        if assigned_class == 0:
                            cm.add_fp()
                        else:
                            cm.add_tp()

                    # set predicted class
                    gene_prediction.set_predicted_class(predicted_class=predicted_class, log_p=log_p)

                    prediction_results[gnid] = gene_prediction

                # set confusion matrix for each fold
                dataset.set_confusion_matrix(cm=cm, fold_idx=fold_idx)

        # TEST
        print ('prediction results length: ', len(prediction_results))

        #return sorted(prediction_results)
        return prediction_results

    def get_fold_size(self):
        return self.fold_size()

