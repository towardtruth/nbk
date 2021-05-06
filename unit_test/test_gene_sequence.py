import unittest
from models import GeneSequence
from utils.DBManager import Pgsql
import settings


class TestGeneSequence(unittest.TestCase):

    # sample sequences
    seq1 = '''CTACGTTACCCATAATATGTGGAAGTATCGACGTAATTTCATAGAGTCATTCGATCTGAATGCTACATGAAGAACATAAGCCAGATGACGGAACGCGGAGACCTAGGATGTAGAAGATCATAACATGAGCGATTCGGCGGATTTGGATTCCTTTTCTATATATCCACTCATGTGGTACTTCATCATACGATTCATATAAGATCCATCTGTCTAGAGATCGTCATATACATCTAGAAAGCCGTATGCTTTGGAAGAAGCTTGTACAGTTTGGGAAGGGGTTTTTTGAGAAAAAAGAAGAATCTACTTCAACCGATATGCCCTTAGGCACGGCCATACATAACATAGAAATCACACGTGGAAGGGGTGGGCAATTAGCTAGAGCAGCAGGTGCTGTAGCGAAACTGATTGCAAAAGAGGGTAAATTGGCCACTTTAAGATTACCATCTGGGGAGGTCCGTTTGGTATCCCAAAACTGCTTAGCAACAGTCGGACAAGTGGGTAATGTTGGGGTGAACCAAAAAAGTTTGGGTAGAGCCGGATCTAAGTGTTGGCTAGGTAAACGCCCCGTAGTAAGAGGGGTAGTTATGAACCCTGTGGACCACCCCCATGGGGGTGGTGAAGGGAAAGCCCCCATTGGTAGAAAAAAACCCACAACCCCTTGGGGTTATCCTGCGCTTGGAAGAAGAACTAGGAAAAGGAAAAAATATAGTGATAGTTTTATTCTTCGTCGCCGTAAGTAAATACGTAACTAGGAATATGGAAAATTGCATTTTTGGAATTTGCAATAATGCGATGGGCGAACGACGGGAATTGAACCCGCGCATGGTGGATTCACAATCCACTGCCTTGATCCACTTGGCTACATCCGCCCCTTATCCAGCTAAAGGATTTTCTCTTTTTTCCATTCATTATTATTCTATTTATTCTGACCTCCATACCTCGATCGAGATAGTGGACATAGGATGCCACTCTTTAAAATGAAAAAAAGGAGTAATCAGCTGTGACACGAAAAAAAACGAATCCTTTTGTAGCTCGTCATTTATTGGCAAAAATCGAAAAGGTCAATATGAAGGAGGAGAAAGAAATAATAGTAACGTGGTCCCGGGCATCTAGCATTCTACCCGCAATGGTTGGCCATACAATCGCGATTCATAATGGAAAAGAACATATACCTATTTACATAACAAATCCTATGGTAGGTCGCAAATTGGGGGAATTCGTGCCTACTCGGCATTTCACGAGTTATGAAAGTACAAGAAAGGATACTAAATCTCGTCGTTAATTGAATTCAGAA'''

    def test_gene_sequence(self):
        conn = Pgsql.Common.connect(settings.conn_string_test)
        #gnid = 58737
        gnid = 76
        seq_type = 'm1'
        gs_pep = GeneSequence(gnid=gnid, seq_type=seq_type, is_max_seq_len=True, conn=conn)
        print(gs_pep.get_seq_str())
        kf = gs_pep.get_kmer_freq(k=7)
        kf.print(sort_type=4, limit=10)


if __name__ == '__main__':
    unittest.main()
