#
# Settings
#

import platform

################################################################################
#
# CONSTANTS
#
################################################################################

#
# Version
#
DEV_VERSION = "1.23.3.25526"


#
# Gene assignment (POS/NEG)
#
gene_id_pos = list()
gene_id_neg = list()

#
# database
#

# PostgreSQL
conn_string = "host='kyoung-z620.usda.iastate.edu' dbname='gene_prediction_p1' user='kyoung' password='**********'"
# use for TEST DB
conn_string_test = "host='kyoung-z620.usda.iastate.edu' dbname='gene_prediction_t1' user='kyoung' password='**********'"


#
# TEST MODE
#
# TEST MODE - Random assigned genes
TEST_MODE_FALSE = 0
TEST_MODE_TRUE = 1
TEST_MODE_NORM = 1
TEST_MODE_RND = 2
TEST_MODE_KMER_FREQ = 11

#
# Validation
#

# validation mode
VALIDATION_MODE_NORMAL = 1
VALIDATION_MODE_SMALL_GENE = 2

# Negative class Building Mode
NEG_CLASS_MODE_NOT_P = 1    # all other genes rather than positive class
NEG_CLASS_MODE_RND_S = 2    # balanced negative class (randomly picked)
NEG_CLASS_MODE_RND_M = 3    # multiple negative classes
NEG_CLASS_MODE_FROM_FS = 4  # set negative class from feature_set table

# Gene loading mode
GN_LD_MODE_PL = 1       # pre-load
GN_LD_MODE_DL = 2       # dynamic load

#
# Sorting types
#
SORT_BY_FREQ_ASC = 1
SORT_BY_FREQ_DESC = 2
SORT_BY_KMER_ASC = 3
SORT_BY_KMER_DESC = 4
SORT_BY_NOSORT = 5

#
# Features (classes, feature name)
#
FN_PA_001 = {'PA-001', 'pa1'}   # PA data   ALL tissues    Top 10%
FN_PA_002 = {'PA-002', 'pa2'}   # PA data   ALL tissues    Top 15%
FN_GE_001 = {'GE-001', 'ge1'}   # GE data   ALL tissues    Top 10%
FN_HGLP_01 = {'HGLP-01'}        # High GE and Low PA, All tissues, g > 50, p < 700
FN_GPCB = {'GPCB'}              # GE and PA Combo, All tissues, two GP combination
FN_GE_N90 = {'GE-N90'}          # GE data   ALL tissues    TOP 10%     using new percentile
FN_PA_N90 = {'PA-N90'}          # PA data   ALL tissues    TOP 10%     using new percentile
FN_GE_N = {'GE-N'}              # GE data   ALL tissues    TOP percentile with Multiple cutoffs     using new percentile
FN_PA_N = {'PA-N'}              # PA data   ALL tissues    TOP percentile with Multiple cutoffs     using new percentile
FN_GE_B = {'GE-B'}              # GE data   ALL tissues    BOTTOM percentile with Multiple cutoffs     using new percentile
FN_PA_B = {'PA-B'}              # PA data   ALL tissues    BOTTOM percentile with Multiple cutoffs     using new percentile


# sqlite3
DEV_DATA_DIR = 'C:\Dev\Data\maize'
DB_PATH = DEV_DATA_DIR + '\maize_AGPv4_pep_all.db3'
#SRC_DATA_DIR = '\data'
#RESULT_DIR = '\results'
if platform.system() == 'Windows':
    DEV_DATA_DIR = 'C:\Dev\Data\maize'
    DB_PATH = DEV_DATA_DIR + '\maize_AGPv4_pep_all.db3'
    #SRC_DATA_DIR = '\data'
    #RESULT_DIR = '\results'
elif platform.system() == 'Linux':
    DEV_DATA_DIR = '/home/kyoung/Taks/Dev/Data/maize'
    DB_PATH = DEV_DATA_DIR + '/maize_AGPv4_pep_all.db3'
    #SRC_DATA_DIR = '/data'
    #RESULT_DIR = '/results'
elif platform.system() == 'Darwin':
    DEV_DATA_DIR = '/home/kyoung/Taks/Dev/Data/maize'
    DB_PATH = DEV_DATA_DIR + '/maize_AGPv4_pep_all.db3'

#
# Data files
#

# Jesses' expression and abundance data
exp_abd = 'exp_abd_1.csv'


# Walley full data
# gene
walley_gene = 'data/gene.out.csv'
# protein
walley_prot = 'data/prot.out.csv'


#
# Sequence files
#

# Peptide sequence
seq_pep = 'data/Zea_mays.AGPv4.pep.all.fa'
# DNA sequence
seq_dna = 'data/RefGen_v4_upstream_sequence.txt'
# Promoter sequence (DNA)
#seq_pmt = '../data/B73_V4_Gene_Model_w_500bp_upstream.txt'
seq_pmt = 'data/B73_V4_Gene_Model_w_500bp_upstream.txt'
seq_pmt_1k = 'data/Zea_mays.AGPv4.genes+1kb_upstream.fa'
seq_pmt_5k = 'data/Zea_mays.AGPv4.genes+5kb_upstream.fa'

#
# Result directory path
#
#RESULT_DIR = '../results/'
RESULT_DIR = '../../nbk/5_experiments/'
RESULT_DIR_COMB = RESULT_DIR + '0_comb/'
RESULT_DIR_EXP_1 = RESULT_DIR + '0_EXP-01/'
RESULT_DIR_EXP_2 = RESULT_DIR + '0_EXP-02/'
RESULT_AGG_SUM_DIR = RESULT_DIR + 'summary/'
# Arffs for Weka Experiments (EXP5-8)
ARFF_DIR_EXP5 = RESULT_DIR + '0_3_EXP-5/0_arffs/'
ARFF_DIR_EXP6 = RESULT_DIR + '0_3_EXP-6/0_arffs/'
ARFF_DIR_EXP7 = RESULT_DIR + '0_3_EXP-7/0_arffs/'
ARFF_DIR_EXP8 = RESULT_DIR + '0_3_EXP-8/0_arffs/'
# RA EXPs
NBK2_RESULT_DIR = '../../nbk2/5_experiments/'

#
# ReExp for find files
#

# feature vector file
RE_FEATURE_VECTOR_START = '^feature_vector'
RE_FEATURE_VECTOR_DETAIL_INFO = '_[gpb]t[0-9][0-9][-].*[_][0-9]'
RE_VERSION = '[0-9].*[.].*[.].*[.].*[0-9]'
RE_FEATURE_VECTOR_VERSION = '_' + RE_VERSION
# prediction result summary file
RE_SUMMARY_START = '^prediction_summary'
RE_FS_DETAIL_INFO = '^[gpb]t[0-9][0-9][-].*'        # Feature Detail Information
RE_SUMMARY_DETAIL_INFO = RE_FS_DETAIL_INFO
# prediction details file
RE_DETAILS_START = '^predictin_details'
# arff file
RE_ARFF_START = '^maize_gp'
# k_size
RE_KSIZE = ''

#
# FileList
#

# list mode: file only, dir only, and all
FLS_ALL = 0
FLS_FILE_ONLY = 1
FLS_DIR_ONLY = 2
# file/dir type
FD_DIR = 1
FD_FILE = 2
# file extensions
FLS_EXT_ALL = ".*"
FLS_EXT_CSV = ".csv"
FLS_EXT_ARFF = ".arff"
FLS_EXT_TXT = ".txt"


#
# Features
#
#FEATURE_SIZE = 26
FEATURE_SIZE = 23


#
# K-mer
#
KMER_MAX_WINDOW_SIZE = 20


#
# Reduced Alphabet mode
#
RA_MODE = False
RAID = -1
NEED_NEW_GENE_DATA = True
RA_AB_SIZE = 0
RA_RP_NUM = 0
ra_mapping_dict = dict()

