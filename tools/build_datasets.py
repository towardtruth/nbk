'''
    Build Data Set in Database
    file: build_datasets.py
    Author: Kyoung Tak Cho
    Created: Mon Jun 18 15:29:52 CDT 2018
    Updated: Fri Jun 21 13:28:03 CDT 2019
    Updated: Wed Jul  3 03:14:19 CDT 2019
    Updated: Saturday, June 20, 2020 10:10:18 AM
'''

import settings
import argparse
import sys
import utils.FileManager as fManager
import utils.DBManager as dbm

choice_yes = {'yes', 'y', 'Yes', 'Y'}
choice_no = {'no', 'n', 'No', 'N'}

parser = argparse.ArgumentParser()
parser.add_argument('--dataset', nargs='+', default=None,
                    help='''name of data set(s) to build.
                            Possible datasets: walley_data or wd
                            | pep_seq_data or pep_seq or pep
                            | dna_seq_data or dna_seq or dna
                            | seq_kmers or kmers or km
                            | promoter_seq_data or promoter or pm''')
                            #| gene_ids or gids or gi''')
parser.add_argument('--db_type', default='pgsql', help='type of database to use')
parser.add_argument('--get_sample_code', default=None, help='get sample codes')
parser.add_argument('-r', '--use_real_db', default='N',
                    choices=choice_yes|choice_no,
                    help='''Set using real (production) database if YES
                    else use TEST DB (test_gene_prediction). Possible options: Y | N and default is N.''')
#parser.add_argument('--class_assign', nargs='+', default=None,
#                    help='''class assignment (features table)
#                            Possible options: {}
#                            | {}
#                            | {}'''.format(settings.FN_PA_001,
#                                            settings.FN_PA_002,
#                                            settings.FN_GE_001))


def main(argv):
    """
    Main
    :param argv:
    """
    args = parser.parse_args(argv[1:])
    if len(argv) <= 1:
        parser.parse_args(['--help'])
        return

    dataset_wd = {'walley_data', 'wd'}
    dataset_pep = {'pep_seq_data', 'pep_seq', 'pep'}
    dataset_dna = {'dna_seq_data', 'dna_seq', 'dna'}
    dataset_km = {'seq_kmers', 'kmers', 'km'}
    dataset_pm = {'promoter_seq_data', 'promoter', 'pm'}
    dataset_kmp = {'seq_kmers_promoter', 'kmers_pm', 'kmp'}
    dataset_ra = {'rda', 'ra'}

    if args.use_real_db:
        if set(args.use_real_db) & choice_yes:
            print('USING TEST DB: NO (USEING REAL/PRODUCTION DB)')
        else:
            settings.conn_string = settings.conn_string_test
            print('USING TEST DB: YES')

    if args.dataset:

        # Build Walley data
        if set(args.dataset) & dataset_wd:
            print('build dataset: walley_data')
            build_walley_data()

        # Build peptide sequence data
        if set(args.dataset) & dataset_pep:
            print('build dataset: pep_seq_data')
            build_pep_seq_data()

        # Build dna sequence data
        if set(args.dataset) & dataset_dna:
            print('build dataset: dna_seq_data')
            build_dna_seq_data()

        # Build kmers for peptide sequence
        if set(args.dataset) & dataset_km:
            print('build dataset: seq_kmers')
            build_kmers()

        # Build kmers for promoter sequence
        if set(args.dataset) & dataset_kmp:
            print('build dataset: seq_kmers_promoter')
            build_kmers(is_promoter=True)

        # import promoter sequence data (DNA sequence)
        if set(args.dataset) & dataset_pm:
            print('import promoter sequence data')
            #build_pmt_seq_data()
            build_pmt_seq_data(bp_size='1k')
            #build_pmt_seq_data(bp_size='5k')

        # build reduced alphabet mapping table
        if set(args.dataset) & dataset_ra:
            build_ra_mapping_table()

    # Get sample codes
    if args.get_sample_code:
        print('Get sample codes')
        print(get_sample_code())


def build_kmers(is_promoter=False, bp_size=None):
    dbManager = dbm.DBManagerPG(settings.conn_string)
    if is_promoter:
        for k in range(7, 0, -1):
            print('{}-mer creating...'.format(k))
            if bp_size is None or bp_size == '1k':
                dbManager.build_promoter_kmers(k=k)
            elif bp_size == '5k':
                dbManager.build_promoter_kmers(k=k, bp_size='5k')
            else:
                dbManager.build_promoter_kmers(k=k)
    else:
        #dbManager.build_kmers(3)
        #dbManager.build_kmers(2)
        dbManager.build_kmers(1)

def build_kmers_pm():
    dbManager = dbm.DBManagerPG(settings.conn_string)
    dbManager.build_kmers(3)

def build_pep_seq_data():
    # Insert into database
    dbManager = dbm.DBManagerPG(settings.conn_string)
    dbManager.import_pep_seq_data(settings.seq_pep)


def build_dna_seq_data():
    # Insert into database
    dbManager = dbm.DBManagerPG(settings.conn_string)
    dbManager.import_dna_seq_data()

# import promoter data
def build_pmt_seq_data(bp_size=None):
    # Insert into database
    dbManager = dbm.DBManagerPG(settings.conn_string)
    if bp_size is None:
        dbManager.import_pmt_seq_data()
    elif bp_size == '1k':
        dbManager.import_pmt_seq_data(seq_type='m1', pmt_fasta_file_path=settings.seq_pmt_1k)
    elif bp_size == '5k':
        dbManager.import_pmt_seq_data(seq_type='m5', pmt_fasta_file_path=settings.seq_pmt_5k)

def build_walley_data():
    # Read data from csv and build data set into list
    fm = fManager.FileManager()
    ld_prot = fm.csv2list(settings.walley_prot)
    ld_gene = fm.csv2list(settings.walley_gene)
    print (len(ld_prot), len(ld_prot[0]))
    print (len(ld_gene), len(ld_gene[0]))

    # Insert into database
    dbManager = dbm.DBManagerPG(settings.conn_string)
    dbManager.build_walley_data(ld_prot, 'p')
    dbManager.build_walley_data(ld_gene, 'g')


#
# Build Reduced Alphabet Mapping Table
#
def build_ra_mapping_table():
    dbManager = dbm.DBManagerPG(settings.conn_string)
    dbManager.build_ra_mapping_table()


def get_sample_code():
    dbManager = dbm.DBManagerPG(settings.conn_string)
    sample_code = dbManager.get_sample_code()
    return sample_code


if __name__ == '__main__':
    main(sys.argv)





