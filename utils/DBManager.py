'''
    Database Manager
    Author: Kyoung Tak Cho
    Created: Mon May 28 23:43:18 CDT 2018
    Updated: Thu Jun 20 13:40:10 CDT 2019
'''


import psycopg2
import sqlite3
from Bio import SeqIO
import settings
import Kmer as km
import sqls
from ReducedAlphabet import ReducedAlphabet

class Pgsql(object):
    class Common(object):
        @staticmethod
        def connect(conn_str=None):
            if conn_str is None:
                conn_str = settings.conn_string
            conn = psycopg2.connect(conn_str)
            #if conn:
            #    print('connected to {} database.'.format(conn_str))
            #else:
            #    print('connection error.')
            return conn

        @staticmethod
        def select_data(sql, pars=None, cur=None, to_list=False):
            use_local_cur = False
            if cur is None or cur.closed:
                use_local_cur = True
            if use_local_cur:
                conn = psycopg2.connect(settings.conn_string)
                cur = conn.cursor()
            # execute
            try:
                if pars is None:
                    cur.execute(sql)
                else:
                    cur.execute(sql % pars)
            except psycopg2.InterfaceError as e:
                print(cur.status)
                print(e)

            # for debugging
            #print('Query: {}'.format(cur.query))

            res = cur.fetchall()
            if to_list:
                res = Pgsql.Common.to_list(res)
            # close
            if use_local_cur:
                cur.close()
                conn.close()
            return res

        @staticmethod
        def insert_data(sql, par_list=None, conn=None, debug=False):
            use_local_conn = False
            if conn is None:
                conn = psycopg2.connect(settings.conn_string)
                use_local_conn = True
            # cursor (local)
            cur = conn.cursor()
            ret = None

            try:
                if debug:
                    print(sql, par_list)
                # execute
                cur.execute(sql, par_list)

                if sql.find('RETURNING') >= 0:
                    ret = cur.fetchone()[0]
                conn.commit()
            except psycopg2.DatabaseError as e:
                print('### Database ERROR ###\n  @Pgsql.Common.insert_data\n\nError Message:')
                print(e)
                print('### END OF ERROR ###')
            finally:
                cur.close()
                if use_local_conn:
                    conn.close()
                return ret

        @staticmethod
        def to_list(src):  # TODO: building  - converting results from db query to list. - to Common library
            res = list()
            for item in src:
                res.append(item[0])
            return sorted(res)


class DBManagerPG(object):

    def __init__(self, conn_string=None):
        self.conn_string = conn_string
        self.conn = None
        self.cur = None

    def __del__(self):
        self.cur.close()
        self.conn.close()

    def connect(self, conn_string=None):
        if conn_string is None:
            conn_string = self.conn_string
        self.conn = psycopg2.connect(conn_string)
        self.cur = self.conn.cursor()

    def connect_new(self, conn_string):
        conn = psycopg2.connect(conn_string)
        return conn

    def disconnect(self):
        self.cur.close()
        self.conn.close()

    def new_cur(self):
        cur = self.conn.cursor()
        return cur

    def get_sample_code(self):
        self.connect()

        self.cur.execute('select * from sample_code')
        rows = []
        for row in self.cur.fetchall():
            cols = []
            for col in row:
                cols.append(col)
            rows.append(cols)

        self.disconnect()

        return rows

    def get_missing_gnids_in_promoter(self):
        sql = sqls.get_missing_gnids_in_promoter
        res = Pgsql.Common.select_data(sql=sql)
        return Pgsql.Common.to_list(res)            # return as a list() type

    def build_walley_data(self, ls_data, vtype):
        self.connect()

        cnt = 1
        for rows in ls_data[1:]:
            gene_id = rows[0]
            self.cur.execute(sqls.build_gene_ids % gene_id)
            for idx, col in enumerate(rows[1:]):
                if col == 'NA':
                    col = 0.
                self.cur.execute(sqls.build_walley_data, (idx + 1, gene_id, vtype, col))
                if cnt % 1000 == 0:
                    print (cnt)
                    cnt += 1

        self.conn.commit()
        self.disconnect()

    def import_pep_seq_data(self, pep_fasta_file_path=settings.seq_pep):
        self.connect()
        items = []
        for seq_record in SeqIO.parse(pep_fasta_file_path, "fasta"):
            gene_id = seq_record.id.split('_')
            pep_seq = str(seq_record.seq)
            pep_seq_len = len(pep_seq)

            # insert gene_id
            self.cur.execute(sqls.build_gene_ids % gene_id[0])

            # build items list - gene_sequence for pep_seq
            items.append([gene_id[0], gene_id[1], pep_seq, pep_seq_len])

        # insert into database
        self.cur.executemany(sqls.build_gene_sequence_pep, items)
        self.conn.commit()
        self.disconnect()

    def import_dna_seq_data(self, dna_fasta_file_path=settings.seq_dna):

        self.connect()
        items = []
        for seq_record in SeqIO.parse(dna_fasta_file_path, "fasta"):
            gene_id = seq_record.id.split('|')
            seq = str(seq_record.seq)
            seq_len = len(seq)

            if seq != 'Sequenceunavailable':
                # insert gene_id
                self.cur.execute(sqls.build_gene_ids % gene_id[0])

                # build items list - gene_sequence for dna_seq
                items.append([gene_id[0], seq, seq_len])

        # insert into database
        self.cur.executemany(sqls.build_gene_sequence_dna, items)
        self.conn.commit()
        self.disconnect()

    def import_pmt_seq_data(self, seq_type='m', pmt_fasta_file_path=settings.seq_pmt):
        self.connect()
        items = []
        cnt_zm = 0
        cnt_en = 0
        cnt_ot = 0
        cnt_en_sub = 0
        cnt_zm_sub = 0
        cnt_na = 0
        for seq_record in SeqIO.parse(pmt_fasta_file_path, "fasta"):
            gene_ids = seq_record.id.split('|')
            seq = str(seq_record.seq)
            seq_len = len(seq)
            gene_id = gene_ids[0]
            gene_id = gene_id.replace('gene:', '')      # if 'gene:' placed in the most lead of each gene_id string, remove it.
            gene_type = gene_id[0:2]
            sub_gene_id = [None, None]

            # get sub_seq info
            if len(gene_ids) >= 2:
                if gene_type == 'Zm':
                    cnt_zm += 1
                    sub_gene_id = gene_ids[1].split('_')
                    if len(sub_gene_id) == 2:
                        cnt_zm_sub += 1
                elif gene_type == 'EN':
                    cnt_en += 1
                    sub_gene_id = gene_ids[1].split('-')
                    if len(sub_gene_id) == 2:
                        cnt_en_sub += 1
                else:
                    cnt_ot += 1

            if seq == 'Sequenceunavailable':
                cnt_na += 1

            # insert gene_id
            self.cur.execute(sqls.build_gene_ids % gene_id)

            # build items list - gene_sequence for dna_seq
            items.append([gene_id, sub_gene_id[1], seq, seq_len, seq_type])

        # test
        print('zm count: {} (sub {}), ENSRNA count: {} (sub {}), other: {}, NA: {}'.format(cnt_zm, cnt_zm_sub,
                                                                                cnt_en, cnt_en_sub,
                                                                                cnt_ot, cnt_na))

        # insert into database
        self.cur.executemany(sqls.build_gene_sequence_pmt, items)
        self.conn.commit()
        self.disconnect()

    def build_kmers(self, k=3,
                    sql_get_seq=sqls.build_kmers_get_seq,
                    sql_put_kmer=sqls.build_kmers_put_kmers):
        self.connect()

        # set db cursor for inserting new kmers
        cur_put = self.new_cur()

        # get all sequence from gene_sequence table
        self.cur.execute(sql_get_seq)
        cnt = 1
        for row in self.cur.fetchall():
            # row[0]: gsid, row[1]: sequence
            freq = km.Kmer.kmer_freq_acc(seq=row[1], window_size=k, freq=None)
            items = []

            for fq in freq:
                items.append((k, row[0], fq, freq[fq]))

            # execution
            cur_put.executemany(sql_put_kmer, items)

            # print message
            #print('k: {}\tgsid: {}\tFreq. Table size: {}'.format(k, row[0], len(freq)))
            #if cnt % 30 == 0:
            if cnt % 100 == 0:
                print('k: {}\tgsid: {}\tFreq. Table size: {}'.format(k, row[0], len(freq)))
                print('\t\t\t\t\t\t<== %dth row done.' % cnt)
            cnt += 1

        cur_put.close()
        self.conn.commit()
        self.disconnect()

    def build_promoter_kmers(self, k=None, bp_size=None):
        if k is None:
            raise ValueError('k is empty.')

        # default settings for sqls
        sql_get_seq = sqls.build_kmers_get_promoter_seq
        sql_put_kmer = sqls.build_kmers_put_kmers_promoter

        if bp_size == '5k':
            sql_get_seq = sqls.build_kmers_get_promoter_seq_5k

        self.build_kmers(k=k, sql_get_seq=sql_get_seq, sql_put_kmer=sql_put_kmer)

    def build_ra_mapping_table(self):
        ra = ReducedAlphabet()
        items = list()

        #for alphabet_size in range(2, 11):
        for alphabet_size in range(11, 21):
            for r in range(1, 11):
                print("Repetation: {}".format(r))
                ra.shuffle()
                ra.partitioning(size=alphabet_size)
                #pprint.pprint(list(ra.group))
                ra.init_mapping_table()
                items.append(('r', alphabet_size, ra.mapping_dict_str()))
                print(ra.mapping_table)

        # mapping
        sql = sqls.build_ra_mapping_table

        self.connect()
        self.cur.executemany(sql, items)
        self.conn.commit()
        self.disconnect()


class DBManager(object):

    #
    # sqlite3
    #
    def __init__(self, db_path):
        self.db_path = db_path
        #self.conn = None
        #self.cur = None

    def connect(self, db_path):
        conn = sqlite3.connect(db_path)
        #self.cur = self.conn.cursor()
        return conn

    def disconnect(self):
        #self.cur.close()
        self.conn.close()

    def new_cursor(self):
        return self.conn.cursor()

    def import_dna_seq(self, dna_fasta_file_path):

        self.connect()
        cur = self.conn.cursor()
        items = []
        for seq_record in SeqIO.parse(dna_fasta_file_path, "fasta"):
            gene_id = seq_record.id.split('|')
            seq = str(seq_record.seq)
            seq_len = len(seq)

            if seq != 'Sequenceunavailable':
                items.append([gene_id[0], seq, seq_len])

        # for test
        '''
        cnt = 1
        for item in items:
            print(item)
            if cnt >= 100:
                break
            cnt += 1
        '''

        # insert into database
        #self.cur.executemany("INSERT INTO gene_sequence (gene_id, seq, seq_len, seq_type) VALUES (?,?,?,'D')", items)
        cur.executemany("INSERT INTO gene_sequence (gene_id, seq, seq_len, seq_type) VALUES (?,?,?,'D')", items)
        self.conn.commit()
        cur.close()

        self.disconnect()


    #
    # Update k-mers for all sequences (maize_pep_seq)
    #
    def build_kmers(k, seq_type):

        SQL = '''
            SELECT gsID, seq
            FROM gene_sequence
            WHERE
              seq_type = ?
            ORDER BY gsID'''

        # sqls
        SQL_GET = ''
        SQL_PUT = ''

        if seq_type == 'P':
            SQL_GET = '''
              SELECT mgID, seq
              FROM maize_pep_seq
              ORDER BY mgID'''
            SQL_PUT = '''
              INSERT INTO kmers
                (k, mgID, kmer, freq)
              VALUES (?,?,?,?)'''
        elif seq_type == 'D':
            SQL_GET = '''
              SELECT gsID, seq
              FROM gene_sequence
              ORDER BY gsID'''
            SQL_PUT = '''
              INSERT INTO seq_kmers
                (k, gsID, kmer, freq)
              VALUES (?,?,?,?)'''
        else:
            SQL_GET = '''
              SELECT mgID, seq
              FROM maize_pep_seq
              ORDER BY mgID'''
            SQL_PUT = '''
              INSERT INTO kmers
                (k, mgID, kmer, freq)
              VALUES (?,?,?,?)'''

        conn = sqlite3.connect(DB_PATH)
        cur = conn.cursor()
        curIns

        cnt = 1
        for row in cur.execute(SQL_GET):
            freq = getNgramFreq(row[1], k)
            update_items = []

            for fq in freq:
                update_items.append((k, row[0], fq, freq[fq]))

            curIns.executemany(SQL_PUT, update_items)
            if cnt % 100 == 0:
                print ('%dth row done.' % cnt)
            cnt += 1

        conn.commit()
        conn.close()

    #
    # Build
    #

#
# main
#


# Build
#dbManager = DBManager('C:\Dev\Data\maize\maize_AGPv4_pep_all.db3')
#dbManager.import_dna_seq("C:\Dev\Data\maize\RefGen_v4_upstream_sequence.txt")






