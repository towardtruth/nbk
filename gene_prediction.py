'''
   Gene prediction
      Author: Kyoung Tak Cho
      Created: Mon Jun 18 11:13:15 CDT 2018
      Updated: Tue Sep 08 19:24:34 CDT 2020
'''

import argparse
import sys
import settings
import sqls
from utils.DBManager import Pgsql
from GeneGroups import Features
from Data import FeatureVector
from models import FeatureInfo, Configuration, Cutoffs, CombConf, RaMap
from utils.Common import DictTool


seq_type_pep = {'pep'}
seq_type_dna = {'dna'}
seq_type_pmt = {'pmt'}
seq_type_rda = {'rda'}
gp_type_g = {'gene', 'g'}
gp_type_p = {'prot', 'p'}
gp_type_b = {'both', 'b'}
db_type_pg = {'pg'}
db_type_sl = {'sl'}
enable_debug = {'y'}
disable_debug = {'n'}
#choice
choice_yes = {'yes', 'y', 'Yes', 'Y'}
choice_no = {'no', 'n', 'No', 'N'}
# build feature group
feature_group_gl = {'gl'}
feature_group_gh = {'gh'}
feature_group_gh10 = {'gh10'}
feature_group_gt = {'gt'}
feature_group_pl = {'pl'}
feature_group_ph = {'ph'}
feature_group_pt = {'pt'}
# build feature vector
feature_vector = {'arff'}
# validation modes
validation_mode_rg = {'reduced_genes', 'rg'}
validation_mode_st = {'small_train_set', 'st'}
gene_load_mode_pl = {'pre_load', 'pl'}
gene_load_mode_dl = {'dynamic_load', 'dl'}


parser = argparse.ArgumentParser()
parser.add_argument('-s', '--seq_type', nargs='+', default='pep',
                   choices=seq_type_pep|seq_type_dna|seq_type_pmt|seq_type_rda,
                    help='''Choose sequence type. ex) -s=pep | dna | pmt | rda''')
parser.add_argument('-p', '--gp_type', default='gene',
                    choices=gp_type_g|gp_type_p|gp_type_b,
                    help='''Choose label data type. ex) gene or prot or both''')
parser.add_argument('-d', '--db_type', default='pgsql',
                    choices=db_type_pg|db_type_sl,
                    help='''Choose database type to use.
                    Possible db: postgresql or pgsql or pg
                    | sqlite3 or sl3 or sl''')
parser.add_argument('-r', '--use_real_db', default='N',
                    choices=choice_yes|choice_no,
                    help='''Set using real (production) database if YES
                    else use TEST DB (test_gene_prediction). Possible options: Y | N and default is N.''')
parser.add_argument('-e', '--enable_debug', default='n',
                    choices=enable_debug|disable_debug,
                    help='''Set debugging mode. Possible options: yes or y | no or n''')
parser.add_argument('-g', '--feature_group', nargs='+', default=None,
                    choices=feature_group_gl|feature_group_gh|feature_group_gh10|feature_group_gt,
                    help='''feature group. Possible parameters:
                    gene_low_exp or gle or gl
                    | gene_high_exp or ghe or gh
                    | gene_high_exp_t10 or ghe10 or gh10
                    | gene_tissues or gt
                    | protein_low_exp or ple or pl
                    | protein_high_exp or phe or ph
                    | protein_tissue or pt''')
parser.add_argument('-c', '--features', nargs='+', default=None,
                    choices=settings.FN_PA_001|settings.FN_PA_002|settings.FN_GE_001
                            |settings.FN_HGLP_01
                            |settings.FN_GPCB
                            |settings.FN_GE_N90|settings.FN_PA_N90
                            |settings.FN_GE_B|settings.FN_PA_B
                            |settings.FN_GE_N|settings.FN_PA_N,
                    help='''class assignment for features. Possible parameters:
                    {}
                    | {}
                    | {}
                    | {}'''.format(settings.FN_PA_001,
                                   settings.FN_PA_002,
                                   settings.FN_GE_001,
                                   settings.FN_HGLP_01,
                                   settings.FN_GPCB,
                                   settings.FN_GE_N90,
                                   settings.FN_PA_N90,
                                   settings.FN_GE_B,
                                   settings.FN_PA_B,
                                   settings.FN_GE_N,
                                   settings.FN_PA_N,))
parser.add_argument('-C', '--percentile', default=None,
                    help='''percnetile''')
parser.add_argument('-m', '--multi_gp', nargs='+', default=None,
                    help='''set GE&PA Combo configurations.''')
parser.add_argument('-u', '--feature_set', type=int, nargs='+',
                    help='''Feature id (fsid) for experiment''')
parser.add_argument('-f', '--feature_vector', nargs='+', default=None,
                    choices=feature_vector,
                    help='''build feature vector. Possible parameters:
                    arff (for Weka)''')
parser.add_argument('-v', '--validation_mode', nargs='+', default=None,
                    choices=validation_mode_rg|validation_mode_st,
                    help='''set validation mode. Possible parameters:
                    reduced_genes_mode or rg
                    | small_train_set or st ''')
parser.add_argument('-t', '--test_mode', default='n',
                    choices=choice_yes|choice_no,
                    help='''enable test_mode''')
parser.add_argument('-i', '--ignore_zero', default='y',
                    choices=choice_yes|choice_no,
                    help='''ignore zero values for data labeling (class assignment)''')
parser.add_argument('-n', '--neg_class_mode', type=int, default=1,
                    help='''Set mode of building negative class set''')
parser.add_argument('-V', '--version', action='version', version=settings.DEV_VERSION,
                    help='''Show current version.''')
parser.add_argument('-G', '--gene_load_mode', default='pre_load',
                    choices=gene_load_mode_pl|gene_load_mode_dl,
                    help='''Gene information loading mode. Possible options:
                    pre_load | dynamic_load''')


#
# main
#
def main(argv):
    # global exp_setting
    exp_setting = Configuration()
    debug_mode = list()
    reduced_mode = False
    test_mode = False
    seq_type = None
    gene_prot = None
    feature_set = list()
    percentile_range = list()

    args = parser.parse_args(argv[1:])
    if len(argv) <= 1:
        parser.parse_args(['--help'])
        return

    # set version
    exp_setting.set_version(settings.DEV_VERSION)
    # Show Version
    version = exp_setting.get_version()
    print('Version:', version.get_version())

    # Get Cutoffs
    cutoffs = Cutoffs()
    cutoffs.query_cutoffs('95, 0, -5')
    exp_setting.set_cutoffs(cutoffs)
    print('Cutoffs data Initialized.')

    # enable debugging mode
    if args.enable_debug:
        if set(args.enable_debug) & enable_debug:
            #debug_mode = [1, 100000]
            debug_mode = [1, 1000]
            reduced_mode = True

    if args.use_real_db:
        if set(args.use_real_db) & choice_yes:
            print('USING TEST DB: NO (USEING REAL/PRODUCTION DB)')
        else:
            settings.conn_string = settings.conn_string_test
            print('USING TEST DB: YES')

    if set(args.test_mode) & choice_yes:
        exp_setting.set_test_mode(True)
        print('TEST MODE: YES')
    else:
        exp_setting.set_test_mode(False)
        print('TEST MODE: NO')

    # ignore zero values
    if set(args.ignore_zero) & choice_yes:
        exp_setting.set_ignore_null(True)
        print('Ignore zero values: YES')
    else:
        exp_setting.set_ignore_null(False)
        print('Ignore zero values: NO')

    # Gene info loading mode
    if args.gene_load_mode:
        if set(args.gene_load_mode) & gene_load_mode_pl:
            exp_setting.set_gene_load_mode(settings.GN_LD_MODE_PL)
            print('Gene loading mode: pre-load')
        elif set(args.gene_load_mode) & gene_load_mode_dl:
            exp_setting.set_gene_load_mode(settings.GN_LD_MODE_DL)
            print('Gene loading mode: dynamic load')

    # sequence type
    if args.seq_type:
        if set(args.seq_type) & seq_type_pep:
            seq_type = 'p'
            print('sequence type: amino acid (peptide)')
        elif set(args.seq_type) & seq_type_dna:
            seq_type = 'd'
            print('sequence type: DNA')
        elif set(args.seq_type) & seq_type_pmt:
            seq_type = 'm1'
            print('sequence type: Promoter data')
            # set missing gnids in promoter data
            exp_setting.set_missing_gnids_in_promoter()
        elif set(args.seq_type) & seq_type_rda:
            seq_type = 'p'
            print('sequence type: Reduced Alphabet')
            settings.RA_MODE = True
    else:   # default
        seq_type = 'p'
        print('sequence type: amino acid (peptide) - Default')
    exp_setting.set_seq_type(seq_type)

    # gp_type
    gp_type = 'g'
    if args.gp_type:
        if set(args.gp_type) & gp_type_g:
            gp_type = 'g'
            print('gp type: g')
        elif set(args.gp_type) & gp_type_p:
            gp_type = 'p'
            print('gp type: p')
        elif set(args.gp_type) & gp_type_b:
            gp_type = 'b'
            print('gp type: b')
        else:
            gp_type = 'g'
            print('gp type: g (default)')
    exp_setting.set_gp_type(gp_type)

    # assign feature groups
    if args.feature_group:
        if set(args.feature_group) & feature_group_gl:
            print('new feature group: gene low expressed')
            Features.gene_low_exp()
        if set(args.feature_group) & feature_group_gh:
            print('new feature group: gene high expressed, top 5%')
            Features.gene_high_exp()
        if set(args.feature_group) & feature_group_gh10:
            print('new feature group: gene high expressed, top 10%')
            Features.gene_high_exp_t10()
        if set(args.feature_group) & feature_group_gt:
            print('new feature group: gene for each tissue, top 10%')
            Features.gene_tissues()

    # feature set
    if args.feature_set:
        feature_set = args.feature_set
        print('feature set: {}'.format(feature_set))

    # set negative class mode
    if args.neg_class_mode:
        neg_class_mode = args.neg_class_mode
        print('NEG_CLASS_MODE:', neg_class_mode)
        if neg_class_mode in (
                settings.NEG_CLASS_MODE_NOT_P,
                settings.NEG_CLASS_MODE_RND_S,
                settings.NEG_CLASS_MODE_RND_M):
            exp_setting.set_neg_class_mode(neg_class_mode)
        else:
            error_mesg = 'NEG_CLASS_MODE:', neg_class_mode, 'is UNKNOWN.'
            raise ValueError(error_mesg)

    # set percentile for new feature set
    if args.percentile:
        percentile_range = args.percentile.split(',')
        percentile_range = [int(x) for x in percentile_range]   # str -> int type
        print('Set percnetile range:', args.percentile)

    # set gp combo configurations
    feature_set_gp_comb = list()
    if args.multi_gp:
        multi_gp_conf = args.multi_gp
        for conf in multi_gp_conf:
            print(conf)
            conf_list = conf.split(':')
            feature_set_gp_comb.append(conf_list)
        print(feature_set_gp_comb)

    # class assignment for features
    if args.features:

        if set(args.features) & (settings.FN_GE_N | settings.FN_GE_B | settings.FN_PA_N | settings.FN_PA_B):
            if len(percentile_range) <= 0:
                raise ValueError('percentile range is empty. Please set percentile range.')

            '''
                It supports adding multiple features at the same time, so it needs to do independently as belows.
            '''
            if set(args.features) & settings.FN_GE_N:
                print('GE_N')
                is_top = True
                gp_type = 'g'
                feature_set_name = next(iter(settings.FN_GE_N))
                for percentile in range(percentile_range[0], percentile_range[1], percentile_range[2]):
                    add_feature_by_percentile(gp_type=gp_type, feature_set_name=feature_set_name,
                                              percentile=percentile, is_top=is_top)
            if set(args.features) & settings.FN_GE_B:
                print('GE_B')
                is_top = False
                gp_type = 'g'
                feature_set_name = next(iter(settings.FN_GE_B))
                for percentile in range(percentile_range[0], percentile_range[1], percentile_range[2]):
                    add_feature_by_percentile(gp_type=gp_type, feature_set_name=feature_set_name,
                                              percentile=percentile, is_top=is_top)
            if set(args.features) & settings.FN_PA_N:
                print('PA_N')
                is_top = True
                gp_type = 'p'
                feature_set_name = next(iter(settings.FN_PA_N))
                for percentile in range(percentile_range[0], percentile_range[1], percentile_range[2]):
                    add_feature_by_percentile(gp_type=gp_type, feature_set_name=feature_set_name,
                                              percentile=percentile, is_top=is_top)
            if set(args.features) & settings.FN_PA_B:
                print('PA_B')
                is_top = False
                gp_type = 'p'
                feature_set_name = next(iter(settings.FN_PA_B))
                for percentile in range(percentile_range[0], percentile_range[1], percentile_range[2]):
                    add_feature_by_percentile(gp_type=gp_type, feature_set_name=feature_set_name,
                                              percentile=percentile, is_top=is_top)

        if set(args.features) & settings.FN_GPCB:
            print('GE&PA Combination data')
            for conf in feature_set_gp_comb:
                add_feature_gp_comb(conf, exp_setting)

    # build feature vector
    if args.feature_vector:
        intervals = [1000]
        if set(args.feature_vector) & feature_vector:
            print('build feature vector')
            fs_set_idx = 0
            #build_feature_vector()
            if reduced_mode:
                #for i in range(1,58938, interval):
                for interval in intervals:
                    for i in range(1,39324, interval):
                    #for i in range(16001,39324, interval):
                        debug_mode = [i, interval]
                        exp_setting.set_debug_mode(debug_mode)
                        build_feature_vector(exp_setting)
            else:
                for k in range(3, 8):
                    exp_setting.set_kmer_size(kmer_size=k)
                    exp_setting.set_genes_info(genes_info=None)
                    for fsid in feature_set:
                        # Version Info
                        print('Version:', settings.DEV_VERSION)

                        if fsid == 0:
                            # set feature info with dummy data for small assigned gene at random
                            fs_info = FeatureInfo(fsid=0,
                                                  fs_name='SM_RND',
                                                  gp_type='g',
                                                  class_size=2)
                            # Set assigned genes limit
                            assigned_genes_limit = [int((x + 23 * fs_set_idx) * 10) for x in range(1, 24)]
                            exp_setting.set_assigned_genes_limit(assigned_genes_limit)
                            fs_set_idx += 1
                        else:
                            # get feature set info from DB
                            res_fs_info = Pgsql.Common.select_data(sqls.get_feature_set, (fsid))
                            fs_info = FeatureInfo(fsid=fsid,
                                                  fs_name=res_fs_info[0][0].strip(),
                                                  gp_type=res_fs_info[0][1].strip(),
                                                  class_size=int(res_fs_info[0][2]))
                        exp_setting.set_fs_info(fs_info)

                        # for test
                        print('### MESSAGE ### fsid: {}, fs_name: {}, gp_type: {}, class_size: {}'.format(
                            exp_setting.get_fsid(),
                            exp_setting.get_fs_name(),
                            exp_setting.get_gp_type(),
                            exp_setting.get_class_size()))

                        debug_mode = [1, 0]
                        exp_setting.set_debug_mode(debug_mode)
                        build_feature_vector(exp_setting)

    # single step classification
    if args.validation_mode:
        if set(args.validation_mode) & validation_mode_rg:
            # reduced gene model
            intervals = [1000, 2000, 3000, 4000, 5000]
            print('validation - reduced genes model mode')

            for interval in intervals:
                for i in range(1,39324, interval):
                    #for i in range(16001,39324, interval):
                    debug_mode = [i, interval]
                    exp_setting.set_debug_mode(debug_mode)
                    build_feature_vector(debug_mode=debug_mode,
                                         gene_prot=gene_prot,
                                         seq_type=seq_type)

#
# Store experiment configurations
#
def store_exp_conf(exp_setting=None):
    gene_dataset = exp_setting.get_gene_dataset()
    gene_positive_class = gene_dataset.get_positive_class()
    gene_negative_class = gene_dataset.get_negative_class()


#
# Add Feature set by percentile
#
def add_feature_by_percentile(gp_type=None, feature_set_name=None, percentile=None, is_top=True):
    if gp_type is None:
        raise ValueError('gp_type is empty.')
    if feature_set_name is None:
        raise ValueError('feature_set_name is empty.')
    if percentile is None:
        raise ValueError('Percentile is empty.')

    if is_top:
        top_percentile = 100 - percentile
        top_bottom = 'TOP'
    else:
        top_percentile = percentile
        top_bottom = 'BOTTOM'
    fs_class_size = 2
    fs_name = '{}{}'.format(feature_set_name, percentile)
    description = 'Gene expression data (GE), All (23) tissues, {} {}% with new cutoffs'.format(top_bottom, top_percentile)

    # add feature set
    fsid = Features.add_new_feature_set(fs_name, gp_type, description, fs_class_size)

    print(fs_name, 'class assignment: gene expression data, all 23 tissues, {} {}%'.format(top_bottom, top_percentile))
    Features.class_assign_by_percentile(fsid=fsid, gp_type=gp_type, percentile=percentile, is_top=is_top)
    print('class assigned of', fs_name, 'feature set successfully.')


#
# Add Feature set for GE & PA Combo
#
def add_feature_gp_comb(conf=None, exp_setting=None):
    feature_set_name = conf[0]
    percentiles = [int(conf[1]), int(conf[2])]
    condition = conf[3]
    is_top1 = None
    is_top2 = None
    percentile1 = percentiles[0]
    percentile2 = percentiles[1]
    hl1 = feature_set_name[0].upper()
    hl2 = feature_set_name[2].upper()
    gp1 = feature_set_name[1].lower()
    gp2 = feature_set_name[3].lower()

    # Top_Bottom [1]
    if hl1 == 'H':
        is_top1 = True
    elif hl1 == 'L':
        is_top1 = False
    else:
        raise ValueError('Invalid value of top/bottom type [1].')
    # Top_Bottom [2]
    if hl2 == 'H':
        is_top2 = True
    elif hl2 == 'L':
        is_top2 = False
    else:
        raise ValueError('Invalid value of top/bottom type [2].')

    # GP type [1]
    if not (gp1 == 'g' or gp1 == 'p'):
        raise ValueError('Invalid value of gp type [1].')
    # GP type [2]
    if not (gp2 == 'g' or gp2 == 'p'):
        raise ValueError('Invalid value of gp type [2].')

    #
    # adding feature set
    #
    fs_gp_type = 'b'
    fs_class_size = 2
    fs_name = '{HL1}{GP1}{PT1}_{HL2}{GP2}{PT2}'.format(HL1=hl1, GP1=gp1.upper(), PT1=percentile1,
                                                       HL2=hl2, GP2=gp2.upper(), PT2=percentile2)
    description = 'GE/PA combo data, all 23 tissues: {}'.format(fs_name)
    # new feature set
    fsid = Features.add_new_feature_set(fs_name, fs_gp_type, description, fs_class_size)
    print('FEATURE SET', fs_name, 'ADDED.')

    #
    # adding features
    #
    print('class assignment: ', description)

    # set Comb configuration info
    comb_conf = CombConf(fsid=fsid, fs_name=fs_name, gp_type=fs_gp_type,
                         hl1=hl1, hl2=hl2, gp1=gp1, gp2=gp2,
                         pt1=float(percentile1/100), pt2=float(percentile2/100),
                         condition=condition)
    Features.class_assign_gp_comb(comb_conf, exp_setting)
    print('class assigned successfully. -', fs_name)


#
# build feature vector
#
def build_feature_vector(exp_setting=None):

    if exp_setting.get_gp_type() is None:
        raise ValueError('gene_prot is empty. Please choose gene_prot type either gene expression data or \
                          protein abundance data.')
    if exp_setting.get_seq_type() is None:
        raise ValueError('seq_type is empty. Please set sequence type: DNA or peptide sequence')

    # set target features
    target_features = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23'
    #target_features = '19,20,21'
    exp_setting.set_target_features(target_features=target_features)

    # set class size
    exp_setting.set_class_size(class_size=2)

    # set fold size for cross-validation
    exp_setting.set_fold_size(fold_size=10)

    # Get mapping table if RA_MODE
    mapping_table = list()
    if settings.RA_MODE:
        # sql = sqls.get_ra_mapping_table
        sql = sqls.get_ra_mapping_table_test_per_each
        mapping_table = Pgsql.Common.select_data(sql)
    else:
        #mapping_table.append([0, "", "", "NON_RA_MODE"])
        mapping_table.append([0, "", "", "{}"])

    RepNum = 1
    for row in mapping_table:
        mapping_dict = DictTool.str2dict(row[3])

        if settings.RAID != row[0]:
            if settings.RA_MODE:
                settings.RAID = row[0]
            settings.NEED_NEW_GENE_DATA = True
        else:
            settings.NEED_NEW_GENE_DATA = False

        settings.RA_AB_SIZE = row[2]
        settings.RA_RP_NUM = RepNum
        settings.ra_mapping_dict = mapping_dict

        # Feature Vector
        f_vector = FeatureVector(exp_setting)

        if exp_setting.get_test_mode():
            f_vector.test_features_dataset()
        else:
            #
            # Prediction for each feature using Cross validation
            #
            # 1. cross validation
            f_vector.cross_validation()

            # 2. write prediction results into a file (.csv)
            f_vector.write_prediction_results()
            f_vector.build_feature_vector()
            f_vector.write_feature_vector()
            f_vector.create_arff()
            f_vector.write_prediction_summary()

            # increase version (build) number
            new_version = settings.DEV_VERSION.split('.')
            build_number = int(new_version[3])
            build_number += 1
            new_version[3] = str(build_number)
            settings.DEV_VERSION = '.'.join(new_version)
            RepNum += 1


#
# Run main program
#
if __name__ == '__main__':
    main(sys.argv)






