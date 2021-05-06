'''
    Data Builder
    Author: Kyoung Tak Cho
    Created: Mon June 3 22:21:08 CDT 2018
    updated: Tue Sep  4 01:23:15 CDT 2018
'''

import settings
from utils.DBManager import Pgsql
from GeneGroups import Features, Genes
from models import ValidationDataset, ConfusionMatrix
from Validation import CrossValidation
from utils.Common import CsvTool as ct
from utils.Common import ListTool
from utils.Common import FileManager, VersionManager
import sqls


class GetData(object):

    @staticmethod
    def wd_all_gnid(gene_prot=None, debug_mode=False):
        '''
        Get all gnids from walley data
        :param gene_prot: string
            'g': gene express data
            'p': protein abundance data
        :return:
            all_gnids: list
        '''
        if gene_prot is None:
            error_mesg = 'gene_prot type is empty.'
            raise ValueError(error_mesg)
        #if seq_type is None:
        #    error_mesg = 'sequence type (seq_type) is empty.'
        #    raise ValueError(error_mesg)

        if debug_mode:
            if debug_mode[1] <= 0:
                limit = 'ALL'
            else:
                limit = debug_mode[1]
            offset = debug_mode[0] - 1

            if gene_prot == 'b':
                sql = sqls.get_wd_all_gnid_debug_both
                pars = (limit, offset)
            else:
                sql = sqls.get_wd_all_gnid_debug
                pars = (gene_prot, limit, offset)
        else:
            if gene_prot == 'b':
                sql = sqls.get_wd_all_gnid_both
                pars = None
            else:
                sql = sqls.get_wd_all_gnid
                pars = (gene_prot)

        if pars is not None:
            print(sql % pars)
            all_gnids = Pgsql.Common.select_data(sql, pars)
        else:
            print(sql)
            all_gnids = Pgsql.Common.select_data(sql)
        return GetData.to_list(all_gnids)


    @staticmethod
    def wd_all_gnid_per_tissue(exp_setting, scid):
        ignore_null = '>'
        if not exp_setting.is_ignore_null():
            ignore_null = '>='

        if exp_setting.get_gp_type() == 'b':
            sql = sqls.get_wd_all_gnid_per_tissue_both
            pars = (scid, ignore_null)
        else:
            sql = sqls.get_wd_all_gnid_per_tissue
            pars = (exp_setting.get_gp_type(), scid, ignore_null)

        print(sql % pars)
        all_gnid_per_tissue = Pgsql.Common.select_data(sql, pars)

        return GetData.to_list(all_gnid_per_tissue)


    @staticmethod
    def to_list(src):  # TODO: building  - converting results from db query to list. - to Common library
        res = list()
        for item in src:
            res.append(item[0])
        return sorted(res)


class FeatureVector(object):

    def __init__(self, exp_setting):
        if exp_setting.get_target_features() is None:
            error_mesg = 'no target features string.'
            raise ValueError(error_mesg)

        self.feature_vector = dict()
        self.wd_all_gnids = None
        self.features_str = exp_setting.get_target_features()
        self.features = None
        self.genes = exp_setting.get_genes_info()
        self.conn = None
        self.gene_prot = exp_setting.get_gp_type()
        self.seq_type = exp_setting.get_seq_type()
        self.debug_mode = exp_setting.get_debug_mode()
        self.test_mode = exp_setting.get_test_mode()
        self.fold_size = exp_setting.get_fold_size()
        self.kmer_size = exp_setting.get_kmer_size()
        self.class_size = exp_setting.get_class_size()
        self.exp_setting = exp_setting

        #
        # initialization
        #

        # open database connection
        self.pg_db_connect()
        # set gnids from walley data where gene data (gene_prot='g')
        if not exp_setting.get_test_mode():
            self._set_wd_all_gnids()
        # set target features list
        self._set_features()
        # set groups
        #self.set_groups()
        # set all genes information of Walley data
        if not exp_setting.get_test_mode():
            self._set_genes_info()

    def __del__(self):
        # close database connection
        if self.conn is not None:
            self.conn.close()

    def pg_db_connect(self):
        self.conn = Pgsql.Common.connect()

    def _set_genes_info(self):
        self.genes = self.exp_setting.get_genes_info()
        if self.genes is None:
            self.genes = Genes(gene_ids=self.wd_all_gnids, seq_type=self.seq_type, k=self.kmer_size)
            self.exp_setting.set_genes_info(self.genes)

    def _set_wd_all_gnids(self):
        '''
        Get all gnids from database and
        Set to the self.wd_all_gnids
        :param gene_prot:
            'g': gene express data
            'p': protein abundance data
        :return:
            None
        '''
        self.wd_all_gnids = list()
        #self.wd_all_gnids = GetData.wd_all_gnid(gene_prot=self.gene_prot,
        #                                        seq_type=self.seq_type)
        self.wd_all_gnids = GetData.wd_all_gnid(gene_prot=self.gene_prot,
                                                debug_mode=self.debug_mode)
        #print(len(self.wd_all_gnids))

    #def set_features(self, target_features=None):
    def _set_features(self):
        '''
        Set features
            self.features: Features
        :return:
            self.features
        '''
        if self.features_str is None:
            error_mesg = "no target_feature string."
            raise ValueError(error_mesg)
        #gnid_min_max = [min(self.wd_all_gnids), max(self.wd_all_gnids)]
        gnid_min_max = self._set_min_max()
        self.features = Features(target_features=self.features_str,
                                 gnid_min_max=gnid_min_max,
                                 test_mode=self.test_mode,
                                 exp_setting=self.exp_setting)

        return self.get_features()

    def _set_min_max(self):
        gnid_min_max = [0, 9999999999]
        if self.wd_all_gnids is not None:
            gnid_min_max = [min(self.wd_all_gnids), max(self.wd_all_gnids)]
        return gnid_min_max

    def get_features(self):
        '''
        Return self.features
        :return:
            self.features: list[] * FEATURE_SIZE
        '''
        return self.features

    def init_feature_dataset(self, fold_size=None):
        #for feature in self.features:
        #    print (feature.name)
        print (fold_size)
        self.features.set_feature_dataset(wd_all_gnids=self.wd_all_gnids,
                                          fold_size=fold_size)

    def _print_feature_dataset(self):
        for feature in self.features:
            sub_total = 0
            for dataset in feature.dataset:
                print('# of gnids for each dataset: ', len(dataset))
                sub_total += len(dataset)
                #for fold in dataset.group:
                for fold_num, fold in enumerate(dataset):
                    print ('# of gnids for each fold: ', len(fold))
                    print (fold[:10])
            total = len(self.wd_all_gnids) - sub_total
            print(total)

    def validation_small_genes(self, class_size=None, gene_size=None, kmer_size=None):
        print('Validation - Small genes model')
        for feature in self.features:
            print('Feature:', feature.name)
            #sg_validation = SmallValidation(genes=self.genes,
            #                                all_gnids=self.wd_all_gnids,
            #                                assigned_gnids=feature.assigned_genes,
            #                                class_size=class_size,
            #                                gene_size=gene_size,
            #                                kmer_size=kmer_size)
            #sg_validation.build_dataset(mode=None)
            #prediction_results = sg_validation.validation()
            #feature.set_prediction_results(prediction_restuls=prediction_results)

            # store confusion matrix in each feature
            #feature.set_confusion_matrix_set(cm_set=self.set_confusion_matrix(validation=sg_validation,
            #                                                                  fold_size=sg_validation.fold_size))

    def cross_validation(self):
        # for each feature, build validation groups with fold size
        print('Cross-Validation')
        for feature in self.features:
            print('Feature Name: {}, k-mer size: {}'.format(feature.name, self.kmer_size))

            # get wd_all_gnids per tissue
            scid = feature.corresp_tissue
            #wd_all_gnids_per_tissue = GetData.wd_all_gnid_per_tissue(self.exp_setting, scid)

            cr_validation = CrossValidation(genes=self.genes,
                                            #all_gnids=wd_all_gnids_per_tissue,
                                            all_gnids=None,
                                            class_size=self.class_size, fold_size=self.fold_size,
                                            kmer_size=self.kmer_size,
                                            exp_setting=self.exp_setting)

            cr_validation.build_datasets(assigned_genes=feature.assigned_genes,
                                         neg_class_mode=self.exp_setting.get_neg_class_mode(),
                                         corresp_tissue=feature.corresp_tissue)

            # Do validation and get prediction results
            prediction_results = cr_validation.validation()
            # store prediction results in each feature
            feature.set_prediction_results(prediction_restuls=prediction_results)

            # store confusion matrix in each feature
            #feature.set_confusion_matrix_set(cm_set=cm_set)
            feature.set_confusion_matrix_set(cm_set=self.set_confusion_matrix(validation=cr_validation,
                                                                              fold_size=self.fold_size))

    def test_features_dataset(self):
        print('TEST features dataset')
        for feature in self.features:
            #print('Feature Name: {}'.format(feature.name))
            #print('\tassigned genes: {}'.format(feature.assigned_genes))
            # get wd_all_gnids per tissue
            scid = feature.corresp_tissue
            wd_all_gnids_per_tissue = GetData.wd_all_gnid_per_tissue(self.exp_setting, scid)

            cr_validation = CrossValidation(genes=self.genes,
                                            #all_gnids=self.wd_all_gnids,
                                            all_gnids=wd_all_gnids_per_tissue,
                                            class_size=self.class_size, fold_size=self.fold_size,
                                            kmer_size=self.kmer_size)

            cr_validation.build_datasets(assigned_genes=feature.assigned_genes,
                                         neg_class_mode=self.exp_setting.get_neg_class_mode(),
                                         corresp_tissue=feature.corresp_tissue)

            cr_validation.test_datasets()

    def set_confusion_matrix(self, validation=None, fold_size=None):
        # add up Confusion Matrix value for all classes
        cm_set = [None] * fold_size
        for fold_idx in range(0, fold_size):
            cm_set[fold_idx] = ConfusionMatrix()

        for class_num, dataset in enumerate(validation.datasets):
            for fold_idx, fold in enumerate(dataset):
                #if class_num == 0:
                #    cm_set.append(ConfusionMatrix())
                cm_set[fold_idx] += dataset.get_confusion_matrix(fold_idx=fold_idx)

        # add up for overall
        cm_all = ConfusionMatrix()
        for cm in cm_set:
            cm_all += cm
        cm_set.append(cm_all)

        # print Overall Confusion Matrix
        print('Prediction Summary - Overall')
        print('accuracy: ', cm_all.get_accuracy())
        print('recall: ', cm_all.get_recall())
        print('precision: ', cm_all.get_precision())
        print('F-measure: ', cm_all.get_f_measure())
        print('FP rate: ', cm_all.get_fp_rate())
        print('PPC rate: ', cm_all.get_ppc_rate())
        print('Confusion matrix: ')
        print(cm_all.get_tp(), '\t', cm_all.get_fn())
        print(cm_all.get_fp(), '\t', cm_all.get_tn(), '\n')

        return cm_set

    def write_genes_info(self):
        print('Write genes information to file.')
        #fd = open('genes_info.pkl', 'wb')
        #pickle.dump(self.genes.genes, fd)
        #fd.close()
        ct.dict2csv(self.genes.genes, 'genes_info.csv')

    def write_prediction_results(self):
        for feature in self.features:
            #print (feature.name)
            prediction_results = feature.prediction_results
            wt_buffer = list()
            is_header_added = False

            for gnid, gp in prediction_results.items():
                # add header
                if not is_header_added:
                    wt_buffer.append(gp.get_header_list())
                    is_header_added = True

                # build prediction details
                wt_buffer.append(gp.to_list())

            # write prediciton results to CSV file
            file_name = '{RD}{VER}/{IT}/{RS}/prediction_details_{FN}_K{KM}_{VER}.csv'.format(
                RD=settings.RESULT_DIR,
                VER=settings.DEV_VERSION,
                IT=self.debug_mode[1],
                RS=self.debug_mode[0],
                KM=self.kmer_size,
                FN=feature.name.strip())
            ListTool.list2csv(wt_buffer, file_name)
            print('{} saved.'.format(file_name))

    def write_prediction_summary(self):
        wt_buffer = list()

        # add header
        wt_buffer.append(ConfusionMatrix.get_header_list())

        # build list of summary for each feature and its folds
        for feature in self.features:
            for fold_num, cm in enumerate(feature.cm_set):
                cm_summary = list()
                if settings.RA_MODE:
                    cm_summary.append(settings.RAID)
                else:
                    cm_summary.append(self.debug_mode[0])
                cm_summary.append(self.kmer_size)
                cm_summary.append(feature.name.strip())
                if fold_num == len(feature.cm_set) - 1:
                    cm_summary.append('overall')
                else:
                    cm_summary.append(fold_num + 1)
                if settings.RA_MODE:
                    cm_summary.append(settings.RA_AB_SIZE)
                cm_summary += cm.to_list()
                wt_buffer.append(cm_summary)

        # write summary
        file_name = '{RD}{VER}/{IT}/prediction_summary_K{KM}_{VER}.csv'.format(
            RD=settings.RESULT_DIR,
            KM=self.kmer_size,
            VER=settings.DEV_VERSION,
            IT=self.debug_mode[1])
        ListTool.list2csv(wt_buffer, file_name, mode='a+')

    def write_feature_vector(self):
        for feature in self.features:
            file_name = '{RD}{VER}/{IT}/{RS}/feature_vector_{FN}_K{KM}_{VER}.csv'.format(
                RD=settings.RESULT_DIR,
                VER=settings.DEV_VERSION,
                KM=self.kmer_size,
                IT=self.debug_mode[1],
                RS=self.debug_mode[0],
                FN=feature.name.strip())
            f = FileManager.file_open(file_name, 'w')

            # set write_buffer
            wt_buffer = list()
            is_header_added = False

            # gnids
            gnids = self.exp_setting.get_gene_dataset_gnids_list(feature_id=feature.corresp_tissue)
            for gnid in gnids:
                vector = self.feature_vector.get(gnid, None)
                if vector is not None:
                    # TODO: add header

                    ## add header
                    #if not is_header_added:
                    #    wt_buffer.append(gp.get_header_list())
                    #    is_header_added = True

                    line = ",".join(str(value) for value in vector)
                    predicted_results = feature.prediction_results.get(gnid)
                    if predicted_results is None:
                        data_label = '?'
                    else:
                        data_label = predicted_results.get_assigned_class()
                    f.write("%s,%s,%s\n" % (gnid, line, data_label))
            f.close()

    def create_arff(self):
        # get feature names
        feature_names = list()
        version = VersionManager(version_str=settings.DEV_VERSION)
        for feature in self.features:
            feature_names.append(feature.name.strip())

        # write arff file for each feature
        for feature in self.features:
            file_name = 'maize_gp_{BD}_K{KM}_{FN}'.format(
                BD=version.build,
                KM=self.kmer_size,
                FN=feature.name.strip())
            file_path = '{RD}{VER}/{IT}/{RS}/{FI}.arff'.format(
                RD=settings.RESULT_DIR,
                VER=settings.DEV_VERSION,
                KM=self.kmer_size,
                IT=self.debug_mode[1],
                RS=self.debug_mode[0],
                FI=file_name)
            f = FileManager.file_open(file_path, 'w')

            # arff header
            #f.write('@relation maize-gp-%s\n' % feature.name.strip())
            f.write('@relation %s\n' % file_name)
            for feature_attr in feature_names:
                f.write('@attribute %s numeric\n' % feature_attr)
            f.write('@attribute class {1,0}\n')
            f.write('@data\n')

            # arff data
            # gnids
            gnids = self.exp_setting.get_gene_dataset_gnids_list(feature_id=feature.corresp_tissue)
            for gnid in gnids:
                vector = self.feature_vector.get(gnid, None)
                if vector is not None:
                    line = ",".join(str(value) for value in vector)
                    predicted_results = feature.prediction_results.get(gnid)
                    if predicted_results is None:
                        data_label = '?'
                    else:
                        data_label = predicted_results.get_assigned_class()
                    f.write("%s,%s\n" % (line, data_label))
            f.close()

    #
    # Build feature vectors
    #
    def build_feature_vector(self):
        print('Build Feature Vector')

        #gene_dataset = self.exp_setting.get_gene_dataset()
        all_gnids_list = self.exp_setting.get_gene_dataset_all_gnids_list()

        #for gnid in gene_dataset.get_all_gnids():
        for gnid in all_gnids_list:
        #for gnid in self.wd_all_gnids:
            gnid_vector = list()
            for feature in self.features:
                #if feature.prediction_results is not None:
                if feature.prediction_results:
                    predicted_results = feature.prediction_results.get(gnid)
                    if predicted_results is None:
                        gnid_vector.append('?')
                    else:
                        #gnid_vector.append(feature.prediction_results[gnid].get_predicted_class())
                        gnid_vector.append(predicted_results.get_predicted_class())
            self.feature_vector[gnid] = gnid_vector


class DataSet(object):
    def __init__(self, dsid=None, tissue=None, exp_type=None, exp_level=None, name=None,
                 pos_genes=None, neg_genes=None):
        self.dsid = None
        self.tissue = None
        self.exp_type = None
        self.exp_level = None
        self.name = None
        self.pos_genes = None       # list(), pos_genes_str() - ListTool.list2str()
        self.neg_genes = None

        if dsid is None:
            raise ValueError("dsid is empty.")
        if tissue is None:
            raise ValueError("tissue is empty.")

        # init
        self._set_dataset(dsid, tissue, exp_type, exp_level, name, pos_genes, neg_genes)

    def _set_dataset(self, dsid, tissue, exp_type, exp_level, name, pos_genes, neg_genes):
        self.dsid = dsid
        self.tissue = tissue
        self.exp_type = exp_type
        self.exp_level = exp_level
        self.name = name
        self.pos_genes = pos_genes
        self.neg_genes = neg_genes

    def get_name(self, type=None):
        name = None
        if type is None:
            name = self.name
        elif type == 1:  # type-level-tissue
            name = "{}-{}-{:02d}".format(self.exp_type, self.exp_level, self.tissue)
        elif type == 2: # type-level-tissue
            name = "{}-{:02d}-{}".format(self.exp_type, self.tissue, self.exp_level)
        elif type == 3: # level-type-tissue
            name = "{}-{}-{:02d}".format(self.exp_level, self.exp_type, self.tissue)
        elif type == 4: # level-type-tissue
            name = "{}-{:02d}-{}".format(self.exp_level, self.tissue, self.exp_type)
        elif type == 5: # level-type-tissue
            name = "{:02d}-{}-{}".format(self.tissue, self.exp_type, self.exp_level)
        elif type == 6: # level-type-tissue
            name = "{:02d}-{}-{}".format(self.tissue, self.exp_level, self.exp_type)

        return name

class DataSetContainer(object):
    def __init__(self, exp_types=None, exp_levels=None, tissues=None):
        self.datasets = list()
        self.exp_types = None
        self.exp_levels = None
        self.tissues = None

        self._set_dataset(exp_types, exp_levels, tissues)

    def __iter__(self):
        for ds in self.datasets:
            yield ds

    def _set_dataset(self, exp_types=None, exp_levels=None, tissues=None):
        if exp_types is None:
            self.exp_types = "'RA','PA'"
        else:
            self.exp_types = exp_types
        if exp_levels is None:
            self.exp_levels = "'T95', 'T90', 'T85', 'T80', 'T75', 'T70'"
        else:
            self.exp_levels = exp_levels
        if tissues is None:
            self.tissues = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23'
        else:
            self.tissues = tissues

        pars = (self.exp_types, self.exp_levels, self.tissues)

        conn_kmd = Pgsql.Common.connect(settings.conn_string_kmd)
        res = Pgsql.Common.select_data(sql=sqls.get_dataset, pars=pars, cur=conn_kmd.cursor())
        for row in res:
            dsid = row[0]
            tissue = row[1]
            e_types = row[2].rstrip()
            e_levels = row[3].rstrip()
            name = row[4].rstrip()
            pos_class = row[5]
            neg_class = row[6]
            self.datasets.append(DataSet(dsid, tissue, e_types, e_levels, name, pos_class, neg_class))

        conn_kmd.close()


class GeneSeq(object):
    def __init__(self, gsid=None, seq=None, gnid=None):
        self.gsid = None
        self.seq = None
        self.gnid = None

        if gsid is None:
            raise ValueError("gsid is empty.")
        if seq is None:
            raise ValueError("seq is empty.")
        if gnid is None:
            raise ValueError("gnid is empty.")

        self.gsid = gsid
        self.seq = seq
        self.gnid = gnid


class GeneSeqContainer(object):
    def __init__(self, gsids_str=None):
        self.gene_seqs = list()
        self.gsids_str = None

        if gsids_str is None or len(gsids_str) == 0:
            raise ValueError("gsids_str is empty.")

        self.gsids_str = gsids_str

        # set gene_seqs
        self._set_gens_seqs()

    def __iter__(self):
        for gene_seq in self.gene_seqs:
            yield gene_seq

    def _set_gens_seqs(self):
        conn_kmd = Pgsql.Common.connect(settings.conn_string_kmd)
        seqs = Pgsql.Common.select_data(sql=sqls.get_seq_by_gsids,
                                        pars=(self.gsids_str), cur=conn_kmd.cursor())
        for seq in seqs:
            gsid = seq[0]
            gnid = seq[1]
            seq = seq[2]
            self.gene_seqs.append(GeneSeq(gsid, seq, gnid))
            #print(gsid)
            #print(gnid)
            #print(seq)

        conn_kmd.close()



