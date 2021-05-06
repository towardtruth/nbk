'''
    Gene Group manager for building dynamic dataset
    file: GeneGroups.py
    Author: Kyoung Tak Cho
    Created: Mon Jul 23 23:25:05 CDT 2018
    Updated: Sun Aug 19 11:36:24 CDT 2018
'''

import settings
import utils.DBManager as dbm
from utils.DBManager import Pgsql
from utils.Common import ListTool
import sqls
from itertools import chain
from models import Feature, Gene, GeneSequence, FreqDict


class FeatureSet(object):
    def __init__(self, name=None, fsid=None, label_data_type=None, description=None,
                 class_size=None, features=None, tissues=None):
        self.name = None
        self.fsid = None
        self.label_data_type = None
        self.description = None
        self.class_size = None
        self.features = None
        self.tissues = None

        # initial set
        self._set(name, fsid, label_data_type, description, class_size, features, tissues)

    def _set(self, name=None, fsid=None, label_data_type=None, description=None,
             class_size=None, features=None, tissues=None):
        if name is not None:
            self.name = name
        if fsid is not None:
            self.fsid = fsid
        if label_data_type is not None:
            self.label_data_type = label_data_type
        if description is not None:
            self.description = description
        if class_size is not None:
            self.class_size = class_size
        if features is not None:
            self.features = features
        if tissues is not None:
            self.tissues = None

    @staticmethod
    def get(fsid=None):
        if fsid is None:
            error_mesg = 'Error: Feature Set ID (fsid) is empty.'
            raise ValueError(error_mesg)
        res = Pgsql.Common.select_data(sqls.get_feature_set, (fsid))
        fs = FeatureSet(name=res[0][0], fsid=fsid, label_data_type=res[0][1], class_size=res[0][2])
        return fs

    @staticmethod
    def tissues(fsid=None):
        if fsid is None:
            error_mesg = 'Error: Feature Set ID (fsid) is empty.'
            raise ValueError(error_mesg)
        tissue_ids = list()
        res = Pgsql.Common.select_data(sqls.get_crsp_tissues_fsid, (fsid))
        tissue_ids = Pgsql.Common.to_list(res)

        return tissue_ids

    def exp_level(self):
        if self.name is None:
            error_mesg = 'Error: Feature Set name is empty.'
            raise ValueError(error_mesg)
        if self.label_data_type is None:
            error_mesg = 'Error: Feature Set label_type is empty.'
            raise ValueError(error_mesg)

        if self.label_data_type == 'b':
            return self.name
        else:
            return self.name[3:]

    def exp_type(self):
        if self.label_data_type is None:
            error_mesg = 'Error: Feature Set label_type is empty.'
            raise ValueError(error_mesg)
        return self.label_data_type


class Features(object):
    """ Features class
            Assign feature groups for Walley dataset

        Attributes:
            features: list[] * FEATURE_SIZE
            target_features: string
                ex) '4,9,16,20': four features
    """

    def __init__(self, target_features=None, gnid_min_max=None, test_mode=False, exp_setting=None):
        if target_features is None:
            error_mesg = 'target features is empty.'
            raise ValueError(error_mesg)

        if gnid_min_max is None:
            error_mesg = 'gnid min/max is empty.'
            raise ValueError(error_mesg)

        self.features = [None] * settings.FEATURE_SIZE
        self.target_features_str = target_features      # selected features number (string type)
        self.target_features = list()                   # selected features number (list type)
        self.gnid_min_max = gnid_min_max
        self.exp_setting = exp_setting

        # init features with targets
        self._init_features()
        self._set_features(test_mode=test_mode)

    def __iter__(self):
        self.n = 0
        return self

    def __next__(self):
        if self.n >= len(self.target_features) or self.n < 0:
            raise StopIteration

        feature = self.features[self.target_features[self.n]]
        self.n += 1

        if feature is not None:
            return feature
        else:
            raise StopIteration

    def __len__(self):
        return len(self.target_features)

    def get_features(self):
        '''
        Return self.features
        :return:
            self.features: list[] * FEATURE_SIZE
        '''
        return self.features

    def get_target_features(self):
        return self.target_features

    def get_feature(self, idx):
        return self.features[idx]

    #def _init_features(self, target_features=None):
    def _init_features(self):
        '''
        Initialize features with target features
            self.target_features: list[] * FEATURE_SIZE
            set True on selected features (fgid)
        :param target_features:
            target feature string with ',' separator
            example: '4,9,16,20'
        :return:
            self.target_features: list[] * FEATURE_SIZE
        '''
        #if target_features is None:
        #    print("Error: no target features string.")
        #    return
        features = str(self.target_features_str).split(',')
        for fid in features:
            fIdx = int(fid) - 1
            self.features[fIdx] = list()
            self.target_features.append(fIdx)

        return self.get_features()

    #
    # Get assigned gnids (Positive class)
    #
    @staticmethod
    def get_gnids_pos(fsid=None, tid=None, cur=None):
        if fsid is None:
            error_mesg = 'Error: Feature set id (fsid) is empty.'
            raise ValueError(error_mesg)
        if tid is None:
            error_mesg = 'Error: Tissue id (tid) is empty.'
            raise ValueError(error_mesg)

        res = Pgsql.Common.select_data(sql=sqls.get_features_info_by_fsid_corresp_tissue,
                                       pars=(fsid, tid), cur=cur)
        #feature_name = res[0][0]
        assigned_genes_str = res[0][1]
        #corresp_tissue = res[0][2]
        assigned_gnids = list(map(int, assigned_genes_str.split(',')))

        return assigned_gnids

    #
    # Get assigned gsids (positive class only)
    #
    @staticmethod
    def get_gsids_pos(gnids=None, seq_type=None, get_max_len=True, cur=None):
        if gnids is None:
            error_mesg = 'Error: gnids is empty.'
            raise ValueError(error_mesg)
        if seq_type is None:
            error_mesg = 'Error: seq_type is empty.'
            raise ValueError(error_mesg)

        gsids = list()

        for gnid in gnids:  # we expect that gnids is list() type
            if get_max_len:
                sql = sqls.get_seq_info_by_gnid_seq_len.format(MM='DESC')   # get max(seq_len)
            else:
                sql = sqls.get_seq_info_by_gnid_seq_len.format(MM='ASC')    # get min(seq_len)
            res = Pgsql.Common.select_data(sql=sql, pars=(gnid, seq_type), cur=cur)
            try:
                gsids.append(res[0][0])
            except IndexError as e:
                print(e)
                print(sql % (gnid, seq_type))
                return None

        return gsids

    #
    # Get gnids by Random
    #
    @staticmethod
    def get_gnids_neg(feature_set_id=None, seq_type=None, label_type=None, tissue_id=None, exclusive_gnids=None, size=None):
        if feature_set_id is None:
            raise ValueError('feature set ID is null.')
        if seq_type is None:
            raise ValueError('sequence type is null.')
        if label_type is None:
            raise ValueError('label type is null')
        if tissue_id is None:
            raise ValueError('tissue ID is null')
        if exclusive_gnids is None:
            exclusive_gnids = '0,0'
        if size is None:
            size = len(exclusive_gnids)

        if seq_type == 'm1':
            # remove missing gnids in promoter
            exclusive_gnids = Features.remove_missing_gnids(exclusive_gnids)

        #pars = (label_type, tissue_id, ListTool.list2str(exclusive_gnids, ','), len(exclusive_gnids))
        pars = (label_type, tissue_id, ListTool.list2str(exclusive_gnids, ','), size)
        random_assigned_gene = Pgsql.Common.select_data(sqls.get_gene_tissues_random, pars)
        random_gnids = ListTool.twoD2oneD(random_assigned_gene)
        return random_gnids

    @staticmethod
    def remove_missing_gnids(self, gnids):
        # remove missing gnids for promoter
        missing_gnids = Features.get_missing_gnids_in_promoter()
        new_sub_dataset = ListTool.sub(gnids, missing_gnids)
        #print('all gnids:', len(all_gnids))
        #print('common items:', len(ListTool.common_items(all_gnids, missing_gnids)))
        #print('removed:', len(new_sub_dataset))
        return sorted(new_sub_dataset)

    @staticmethod
    def get_exclusive_gnids(self, gnids, missing_gnids):
        exclusive_gnids = ListTool.add_rm_dup(gnids, missing_gnids)
        return sorted(exclusive_gnids)

    @staticmethod
    def get_missing_gnids_in_promoter():
        sql = sqls.get_missing_gnids_in_promoter
        res = Pgsql.Common.select_data(sql=sql)
        return Pgsql.Common.to_list(res)

    def _set_features(self, test_mode=False):
        #if self.exp_setting.get_fsid() == 0:
        #    pass

        #assigned_genes_limit = [8294, 2189, 4378,
        #                        1300, 1321,  685, 1001,  958, 1147,  869, 1185,  322, 1451,
        #                        1185,  943,  895, 1062, 1017,  589, 1329, 1389, 1375, 1233,
        #                        1353, 1379, 1232]
        assigned_genes_limit = self.exp_setting.get_assigned_genes_limit()
        #assigned_genes_limit = [int(10 + pow(x, 2.96)) for x in range(1, 24)]
        #assigned_genes_limit = [int(10) for x in range(1, 24)]
        #assigned_genes_limit = [int(x*10) for x in range(1, 24)]
        #assigned_genes_limit = [int((x+23)*10) for x in range(1, 24)]
        #assigned_genes_limit = [int((x+23*2)*10) for x in range(1, 24)]
        #assigned_genes_limit = [int((x+23*3)*10) for x in range(1, 24)]
        #assigned_genes_limit = [int((x+23*4)*10) for x in range(1, 24)]
        #assigned_genes_limit = [int((x+23*5)*10) for x in range(1, 24)]
        #assigned_genes_limit = [int((x+23*6)*10) for x in range(1, 24)]
        #assigned_genes_limit = [int((x+23*7)*10) for x in range(1, 24)]
        #assigned_genes_limit = [int((x+23*8)*10) for x in range(1, 24)]
        #assigned_genes_limit = [int((x+23*9)*10) for x in range(1, 24)]
        #assigned_genes_limit = [int((x+23*10)*10) for x in range(1, 24)]
        #assigned_genes_limit = [int(x*500) for x in range(1, 24)]

        for idx, feature in enumerate(self.features):
            if feature is not None:
                # get gnids
                # TODO: cover for fsid == 0 which is the case of samll assgined gene at random
                if self.exp_setting.get_fsid() == 0:    # small gene mode (test mode)
                    feature_name = 'SAG_RND_t%02d' % (idx + 1)
                    corresp_tissue = idx + 1
                    pars = ('g', corresp_tissue, '0,0', assigned_genes_limit[idx])
                    random_assigned_gene = Pgsql.Common.select_data(sqls.get_gene_tissues_random, pars)
                    assigned_gene = Pgsql.Common.to_list(random_assigned_gene)
                else:
                    res = Pgsql.Common.select_data(sqls.get_features_info_by_fsid_corresp_tissue,
                                                   (self.exp_setting.get_fsid(), idx+1,))
                    feature_name = res[0][0]
                    assigned_genes_str = res[0][1]
                    corresp_tissue = res[0][2]
                    assigned_gene = list(map(int, assigned_genes_str.split(',')))

                #if test_mode:
                #    pars = ('g', corresp_tissue, assigned_genes_str, assigned_genes_limit[idx])
                #    #print('Randomly assigned genes,  Feature name: %s' % feature_name)
                #    #print(sqls.get_gene_tissues_random % pars)
                #    random_assigned_gene = Pgsql.Common.select_data(sqls.get_gene_tissues_random, pars)
                #    assigned_gene = Pgsql.Common.to_list(random_assigned_gene)

                #print('random assigned gene: \n{}'.format(assigned_gene))

                if self.gnid_min_max:
                    assigned_gene = [gnid for gnid in assigned_gene if self.gnid_min_max[0] <= gnid <= self.gnid_min_max[1]]
                #print (len(assigned_gene))
                #print(assigned_gene)
                #print('after applying min_max assigned gene: \n{}'.format(assigned_gene))

                self.features[idx] = Feature(name=feature_name,
                                             class_size=self.exp_setting.get_class_size(),
                                             assigned_genes=assigned_gene,
                                             corresp_tissue=corresp_tissue)

    #def set_feature_dataset(self, wd_all_gnids=None, fold_size=None):
    #
    #    for feature in self:
    #        #print(feature.name)
    #        feature.set_dataset(wd_all_gnids=wd_all_gnids, fold_size=fold_size)


    ###########################################################################
    #
    # Static Methods
    # TODO: remove gnid in walley dataset which has no actual sequence
    #   - Done: with removing items in walley_data table
    #
    ###########################################################################
    @staticmethod
    def _low_gene_exp(sql=None):
        if sql is None:
            sql = sqls.sql_low_exp_gene_ids
        low_gene_exp_gene_ids = dbm.Pgsql.Common.select_data(sql, 'g')
        return low_gene_exp_gene_ids

    @staticmethod
    def _new_feature_group(feature_code=None, class_size=2, sql=None, description=None, corresp_tissue=None, gene_prot='g'):
        # get gene ids
        gene_ids_raw = dbm.Pgsql.Common.select_data(sql)
        gene_ids_raw.sort()    # sort by gnid
        gene_ids = list(chain.from_iterable(gene_ids_raw))

        feature_name = "{GP}{FC}".format(GP=gene_prot, FC=feature_code)
        desc = description + sql
        assigned_genes = ",".join(map(str, gene_ids))

        # insert
        dbm.Pgsql.Common.insert_data(sqls.new_feature_group,
                                     (feature_name, class_size, desc, assigned_genes, corresp_tissue))

    @staticmethod
    def _new_features(fsid=None, feature_code=None, sql=None, description=None, corresp_tissue=None, gene_prot='g'):
        # get gene ids
        gene_ids_raw = dbm.Pgsql.Common.select_data(sql)
        gene_ids_raw.sort()    # sort by gnid
        gene_ids = list(chain.from_iterable(gene_ids_raw))
        feature_size = len(gene_ids)

        feature_name = "{GP}{FC}".format(GP=gene_prot, FC=feature_code)
        desc = description + sql
        assigned_genes = ",".join(map(str, gene_ids))

        # insert
        dbm.Pgsql.Common.insert_data(sqls.new_features,
                                     (fsid, feature_name.strip(), desc, assigned_genes, corresp_tissue, feature_size))

    @staticmethod
    def gene_low_exp():     # TODO: need to think about value_type option
        # set arguments
        gene_prot = 'g'
        feature_code = 'l'
        class_size = 2
        sql = sqls.get_low_exp_gene_ids % gene_prot
        description = "name: low expressed genes,\n\nsql: "
        corresp_tissue = None
        # new class assignment
        Features._new_feature_group(feature_code=feature_code, class_size=class_size,
                                        sql=sql, description=description,
                                        corresp_tissue=corresp_tissue, gene_prot=gene_prot)

    @staticmethod
    def gene_high_exp():
        # set arguments
        gene_prot = 'g'
        feature_code = 'h5'
        class_size = 2
        sql = sqls.get_high_exp_gene_ids % (gene_prot, gene_prot)
        description = "name: high expressed genes (top 5%),\n\nsql: "
        corresp_tissue = None
        # new class assignment
        Features._new_feature_group(feature_code=feature_code, class_size=class_size,
                                        sql=sql, description=description,
                                        corresp_tissue=corresp_tissue, gene_prot=gene_prot)

    @staticmethod
    def gene_high_exp_t10():
        # set arguments
        gene_prot = 'g'
        feature_code = 'h10'
        class_size = 2
        sql = sqls.get_high_exp_gene_ids_t10 % (gene_prot, gene_prot)
        description = "name: high expressed genes (top 10%),\n\nsql: "
        corresp_tissue = None
        # new class assignment
        Features._new_feature_group(feature_code=feature_code, class_size=class_size,
                                        sql=sql, description=description,
                                        corresp_tissue=corresp_tissue, gene_prot=gene_prot)

    @staticmethod
    def gene_tissues():
        # set arguments
        gene_prot = 'g'
        class_size = 2
        cutoff = Features.get_cutoff(gene_prot=gene_prot)

        for tissue in range(1, 24):
            feature_code = 't%s' % tissue
            sql = sqls.get_gene_tissues % (gene_prot, tissue, cutoff)
            description = "name: genes by tissue# {TS} (top 10%),\n\nsql: ".format(TS=tissue)
            corresp_tissue = tissue
            # new class assignment
            Features._new_feature_group(feature_code=feature_code, class_size=class_size,
                                            sql=sql, description=description,
                                            corresp_tissue=corresp_tissue, gene_prot=gene_prot)

    @staticmethod
    def prot_tissues(gene_prot=None, class_size=2, fsid=None):
        pass


    @staticmethod
    def class_assign_tissues(gene_prot=None, fsid=None):
        if fsid is None:
            raise ValueError('fsid is empty.')
        #gene_prot = 'p'
        #class_size = 2
        cutoff = Features.get_cutoff(gene_prot=gene_prot)

        for tissue in range(1, 24):
            feature_code = 't%s' % tissue
            sql = sqls.get_gene_tissues % (gene_prot, tissue, cutoff)
            description = "name: protein abundance by tissue# {TS} (top 10%),\n\nsql: ".format(TS=tissue)
            corresp_tissue = tissue
            # new class assignment
            Features._new_features(fsid=fsid, feature_code=feature_code, sql=sql, description=description,
                                        corresp_tissue=corresp_tissue, gene_prot=gene_prot)

    @staticmethod
    def class_assign_hglp(gp_type=None, fsid=None):
        if fsid is None:
            raise ValueError('fsid is empty.')

        for tissue in range(1, 24):
            feature_code = 't%s' % tissue
            sql = sqls.get_gene_hglp.format(SCID=tissue, GCUT=50, PCUT=700) # % (gene_prot, tissue, cutoff)
            description = "name: Combination of high gene expression and low protein abundance by tissue# {TS},\n\nsql: ".format(TS=tissue)
            corresp_tissue = tissue
            # new class assignment
            Features._new_features(fsid=fsid, feature_code=feature_code, sql=sql, description=description,
                                   corresp_tissue=corresp_tissue, gene_prot=gp_type)

    @staticmethod
    def class_assign_by_percentile(fsid=None, gp_type=None, percentile=None, is_top=True):
        if fsid is None:
            raise ValueError('fsid is empty.')
        if gp_type is None:
            raise ValueError('gp_type is empty.')
        if percentile is None or percentile is '':
            raise ValueError('percentile value is empty.')

        for tissue in range(1, 24):
            if is_top:
                top_bottom = 'N'
                inequality = '>='
            else:
                top_bottom = 'B'
                inequality = '<='
            feature_code = 't%02d-%s%02d' % (tissue, top_bottom, percentile)
            percentile_decimal = '0.%02d' % (percentile)
            sql = sqls.get_gene_by_cutoffs.format(PCNT=percentile_decimal, GP=gp_type, TSID=tissue, INEQ=inequality)
            #print(sql)
            description = "name: gene expression with new cutoff (percentile setting) for tissue# {TS},\n\nsql: ".format(TS=tissue)
            corresp_tissue = tissue
            # new class assignment
            Features._new_features(fsid=fsid, feature_code=feature_code, sql=sql, description=description,
                                   corresp_tissue=corresp_tissue, gene_prot=gp_type)

    @staticmethod
    def class_assign_gp_comb(comb_conf=None, exp_setting=None):
        if comb_conf is None:
            raise ValueError('combination configuration info is empty.')

        # CutOffs data
        cutoffs = exp_setting.get_cutoffs()

        for tissue in range(1, 24):
            # get cutoffs info
            cf_key1 = (comb_conf.gp1, tissue, comb_conf.pt1)
            cf_key2 = (comb_conf.gp2, tissue, comb_conf.pt2)
            cutoff1 = cutoffs.cutoff(cf_key1)
            cutoff2 = cutoffs.cutoff(cf_key2)

            if comb_conf.hl1 == 'H':
                inequality1 = '>'
            else:
                inequality1 = '<'

            if comb_conf.hl2 == 'H':
                inequality2 = '>'
            else:
                inequality2 = '<'

            feature_code = 't%02d-%s' % (tissue, comb_conf.fs_name)
            #percentile1 = comb_conf.pt1
            #percentile2 = comb_conf.pt2
            sql = sqls.get_gene_gp_comb.format(GP1=comb_conf.gp1, GP2=comb_conf.gp2,
                                               SCID=tissue, CUT1=cutoff1, CUT2=cutoff2,
                                               INQ1=inequality1, INQ2=inequality2)
            #print(sql)
            description = "description: GE&PA combination for tissue# {TS},\n\nsql: ".format(TS=tissue)
            print(description)
            corresp_tissue = tissue
            # new class assignment
            Features._new_features(fsid=comb_conf.fsid, feature_code=feature_code,
                                   sql=sql, description=description,
                                   corresp_tissue=corresp_tissue, gene_prot=comb_conf.gp_type)

    @staticmethod
    def class_assign_GE_N90(gp_type=None, fsid=None):
        if fsid is None:
            raise ValueError('fsid is empty.')

        for tissue in range(1, 24):
            feature_code = 't%02d-N90' % tissue
            sql = sqls.get_gene_by_cutoffs.format(PCNT='0.90', GP=gp_type, TSID=tissue)
            description = "name: gene expression with new cutoff (percentile setting) for tissue# {TS},\n\nsql: ".format(TS=tissue)
            corresp_tissue = tissue
            # new class assignment
            Features._new_features(fsid=fsid, feature_code=feature_code, sql=sql, description=description,
                                   corresp_tissue=corresp_tissue, gene_prot=gp_type)

    @staticmethod
    def class_assign_PA_N90(gp_type=None, fsid=None):
        if fsid is None:
            raise ValueError('fsid is empty.')

        for tissue in range(1, 24):
            feature_code = 't%02d-N90' % tissue
            sql = sqls.get_gene_by_cutoffs.format(PCNT='0.90', GP=gp_type, TSID=tissue)
            description = "name: Protein abundance with new cutoff (percentile setting) for tissue# {TS},\n\nsql: ".format(TS=tissue)
            corresp_tissue = tissue
            # new class assignment
            Features._new_features(fsid=fsid, feature_code=feature_code, sql=sql, description=description,
                                   corresp_tissue=corresp_tissue, gene_prot=gp_type)

    @staticmethod
    def add_new_feature_set(fs_name, gp_type, description, class_size):
        # insert feature set
        fsid = dbm.Pgsql.Common.insert_data(sqls.new_feature_set,
                                            (fs_name, gp_type, description, class_size))
        return fsid

    @staticmethod
    def get_cutoff(gene_prot):
        cutoff = 0
        if gene_prot == 'g':
            cutoff = 114.26
        elif gene_prot == 'p':
            cutoff = 1471.77

        return cutoff


#
# Genes class - all genes from walley data
#
class Genes(object):
    def __init__(self, gene_ids=None, seq_type='p', k=3):
        self.genes = dict()
        self.seq_type = seq_type
        self.k = k

        # initialization
        self._set_genes(gene_ids)

    def __iter__(self):
        for gene in self.genes:
            yield gene

    def __len__(self):
        return len(self.genes)

    def _set_genes(self, gene_ids):
        if gene_ids is None:
            return
        conn = Pgsql.Common.connect()
        cnt = 0
        for gnid in gene_ids:
            # get sequence information
            #gene_seq_pep_min = GeneSequence(gnid=gnid, seq_type=self.seq_type, is_max_seq_len=False, conn=conn)
            gene_seq_pep_min = None
            gene_seq_pep_max = GeneSequence(gnid=gnid, seq_type=self.seq_type, is_max_seq_len=True, conn=conn, k=self.k)
            gene_seq_dna_max = None

            # set Gene instance for each gnid
            gene = Gene(gnid=gnid,
                            pep_seq_min=gene_seq_pep_min,
                            pep_seq_max=gene_seq_pep_max,
                            dna_seq=gene_seq_dna_max)

            self.genes[gnid] = gene
            cnt += 1

            # for test
            if cnt % 100 == 0:
                print ('%s gene has been added.' % gnid)

        conn.close()   # for single adding

    def get_gene(self, gnid=None):
        if gnid is None:
            error_mesg = 'gnid is empty.'
            raise ValueError(error_mesg)

        gene = self.genes.get(gnid, 0)
        if gene == 0:
            self._set_genes([gnid])
            gene = self.genes[gnid]

        return gene


#
# Frequency Dictionary Set
#   : for managing kmer frequency dictionaries for each class,
#
class FreqDictSet(object):

    def __init__(self):
        self.fd_fold = list()       # FreqDict * KMER_MAX_WINDOW_SIZE * fold_size
        #self.kf_fold_k2 = list()
        self.fd_class = list()      # FreqDict * KMER_MAX_WINDOW_SIZE

        #
        # Initialization
        #

        # create fd_fold with # of KMER_MAX_WINDOW_SIZE
        self._init_fd_fold()


    def _init_fd_fold(self):
        for k in range(0, settings.KMER_MAX_WINDOW_SIZE):
            self.fd_fold.append(list())     # for each k-mer window size
            #self.fd_fold.append(FreqDict(k=k+1))     # for each k-mer window size
            self.fd_class.append(FreqDict(k=k+1))

    def _assert_k(self, k=None):
        if k is None:
            error_mesg = 'k (window size) is empty.'
            raise ValueError(error_mesg)
        if k < 0 or k >= settings.KMER_MAX_WINDOW_SIZE:
            error_mesg = 'invalide k (window size) value.'
            print ('k = ', k)
            raise ValueError(error_mesg)

    def append_kmer_dict(self, k=None, freq_dict=None):
        '''
            Append freq_dict per fold (FreqDict)
        :param k:
        :param freq_dict:
            frequency dictionary for a fold set
        :return:
            None
        '''
        self._assert_k(k)

        self.fd_fold[k-1].append(freq_dict)

    def set_fd_fold(self, fd_fold=None, k=None):
        if fd_fold is None:
            error_mesg = 'kmer frequency fold is empty.'
            raise ValueError(error_mesg)
        self._assert_k(k)

        self.fd_fold[k-1] = fd_fold

    def get_fd_fold(self, k=None, fold_idx=None):
        self._assert_k(k)

        if fold_idx is None:
            error_mesg = 'fold_idx is empty.'
            raise ValueError(error_mesg)

        return self.fd_fold[k-1][fold_idx]

    def get_fd_class(self, k=None):
        self._assert_k(k)

        if len(self.fd_class[k-1]) == 0:
            fd_class = FreqDict(k=k)
            for fd in self.fd_fold[k-1]:
                fd_class += fd
                #print('fd in loop: ', len(fd), fd.total_freq_float())
                #print('fd_class in loop: ', len(fd_class), fd_class.total_freq_float())
            print('fd_class in loop (final): ', len(fd_class), fd_class.total_freq_float())
            self.fd_class[k-1] = fd_class
        return self.fd_class[k-1]

    def get_fd_complement(self, k=None, fold_idx=None):
        self._assert_k(k)

        if fold_idx is None:
            error_mesg = 'fold_idx is empty.'
            raise ValueError(error_mesg)

        #fd_complement = FreqDict(k=k)
        #fd_class = self.get_fd_class(k=k)

        return self.get_fd_class(k=k) - self.get_fd_fold(k=k, fold_idx=fold_idx)



#
# Group Partitioning class
#
class Partition(object):

    def __init__(self, data, group_size):
        self.group_size = group_size
        self.group = [None] * self.group_size
        self.data = data

    def partition(self):
        for idx, value in enumerate(self.data):
            self.group[idx % self.group_size].append(value)
        groups = Group(self.group)
        return self.group


