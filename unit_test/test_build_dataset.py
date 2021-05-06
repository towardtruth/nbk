import unittest
import settings
from Bio import SeqIO


class TestBuildDatasets(unittest.TestCase):
    def test_promoter_data(self):
        pmt_fasta_file_path = settings.seq_pmt_5k
        count = 0

        for seq_record in SeqIO.parse(pmt_fasta_file_path, "fasta"):
            gene_ids = seq_record.id.split('|')
            seq = str(seq_record.seq)
            seq_len = len(seq)
            gene_id = gene_ids[0]
            gene_id = gene_id.replace('gene:', '')      # if 'gene:' placed in the most lead of each gene_id string, remove it.
            gene_type = gene_id[0:2]
            sub_gene_id = [None, None]

            print('before:', gene_id)
            gene_id = gene_id.replace('gene:', '')
            print('after:', gene_id)
            gene_id = gene_id.replace('gene:', '')
            print('again:', gene_id)

            count += 1
            if count >= 100:
                break



if __name__ == '__main__':
    unittest.main()
