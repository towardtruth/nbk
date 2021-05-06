'''
    Build Data Set in Database
    Author: Kyoung Tak Cho
    Created: Mon Jun 18 15:29:52 CDT 2018
    Updated: Mon Jul 23 22:23:20 CDT 2018
'''

import settings
import argparse
import sys
import utils.FileManager as fManager
import utils.DBManager as dbm

parser = argparse.ArgumentParser()
parser.add_argument('--dataset', nargs='+', default=None,
                    help='''name of data set(s) to build.
                            Possible datasets: walley_data_full or walley_data or wd
                            | pep_seq_data or pep_seq or pep
                            | dna_seq_data or dna_seq or dna
                            | seq_kmers or kmers or km''')
parser.add_argument('--db_type', default='psql', help='type of database to use')
parser.add_argument('--get_sample_code', default=None, help='get sample codes')

#
# build
#
def main(argv):
    args = parser.parse_args(argv[1:])
    if len(argv) <= 1:
        parser.parse_args(['--help'])
        return

    # Get sample code
    #print (get_sample_code())

    dataset_wd = {'walley_data_full', 'walley_data', 'wd'}
    dataset_pep = {'pep_seq_data', 'pep_seq', 'pep'}
    dataset_dna = {'dna_seq_data', 'dna_seq', 'dna'}
    dataset_km = {'seq_kmers', 'kmers', 'km'}

    # Build Walley Full data
    if args.dataset and set(args.dataset) & dataset_wd:
        print ('build dataset: walley_data_full')
        build_walley_full_data()

    # Build peptide sequence data
    if args.dataset and set(args.dataset) & dataset_pep:
        print ('build dataset: pep_seq_data')
        build_pep_seq_data()

    # Build dna sequence data
    if args.dataset and set(args.dataset) & dataset_dna:
        print ('build dataset: dna_seq_data')
        build_dna_seq_data()

    # Build kmers for peptide sequence
    if args.dataset and set(args.dataset) & dataset_km:
        print ('build dataset: seq_kmers')
        build_kmers()

    # Get sample codes
    if args.get_sample_code:
        print ('Get sample codes')
        print (get_sample_code())



def build_kmers():
    """ Build kmers dataset for pepide sequence

        Args:
            None

        Return:
            None
    """
    dbManager = dbm.DBManagerPG(settings.conn_string)
    dbManager.build_kmers(3)
    dbManager.build_kmers(2)

def build_pep_seq_data():
    # Insert into database
    dbManager = dbm.DBManagerPG(settings.conn_string)
    dbManager.import_pep_seq_data(settings.seq_pep)

def build_dna_seq_data():
    # Insert into database
    dbManager = dbm.DBManagerPG(settings.conn_string)
    dbManager.import_dna_seq_data()

def build_walley_full_data():
    # Read data from csv and build data set into list
    fm = fManager.FileManager()
    ld_prot = fm.csv2list(settings.walley_prot)
    ld_gene = fm.csv2list(settings.walley_gene)
    print (len(ld_prot), len(ld_prot[0]))
    #print(ls_data[:10])
    print (len(ld_gene), len(ld_gene[0]))

    # Insert into database
    dbManager = dbm.DBManagerPG(settings.conn_string)
    dbManager.build_walley_full_data(ld_prot, 'p')
    dbManager.build_walley_full_data(ld_gene, 'g')

def get_sample_code():
    dbManager = dbm.DBManagerPG(settings.conn_string)
    sample_code = dbManager.get_sample_code()
    return sample_code

if __name__ == '__main__':
    main(sys.argv)





