import unittest
import sqls
from utils.DBManager import Pgsql
import settings
from models import FreqDict


class TestPromoterSeq(unittest.TestCase):
    def test_get_promoter_sequence(self):
        build_kmers_get_promoter_seq = '''
          SELECT gsid, LEFT(seq, 1000)
          FROM gene_sequence
          WHERE
            seq_type = 'm1'
          ORDER BY gsid
          LIMIT 10'''

        conn = Pgsql.Common.connect(settings.conn_string_test)
        cur = conn.cursor()

        #sql = sqls.build_kmers_get_promoter_seq
        sql = build_kmers_get_promoter_seq
        gs_info = Pgsql.Common.select_data(sql=sql, cur=cur)
        cnt = 0
        for gsid, seq in gs_info:
            print(gsid, seq)
            cnt += 1

            if cnt >= 10:
                break
        cur.close()
        conn.close()

    def test_get_kmers(self):
        get_kmer_freq_by_gsid = """
          SELECT kmer, freq
          FROM %s
          WHERE
            k = %s
            AND gsid = %s
          ORDER BY kmer"""

        conn = Pgsql.Common.connect(settings.conn_string_test)
        cur = conn.cursor()

        k = 5
        sql = sqls.build_kmers_get_promoter_seq
        gs_info = Pgsql.Common.select_data(sql=sql, cur=cur)
        cnt = 0
        kmer_freq_sum = FreqDict(k=k)
        for gsid, seq in gs_info:
            pars_kmer = ('kmers_promoter', k, gsid)
            res = Pgsql.Common.select_data(sql=sqls.get_kmer_freq_by_gsid, pars=pars_kmer, cur=cur)
            kmer_freq = FreqDict(k=k)

            for row in res:
                kmer = str(row[0]).upper()
                freq = row[1]
                kmer_freq.add_kmer_freq(kmer=kmer, freq=freq)

            kmer_freq_sum += kmer_freq

            print('gsid: {}, frequency dictionary size: {} and total: {}'.format(gsid, len(kmer_freq), len(kmer_freq_sum)))

            #cnt += 1

            #if cnt >= 1000:
            #    break

        cur.close()
        conn.close()


if __name__ == '__main__':
    unittest.main()
