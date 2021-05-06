'''
    All models
    file: models.py
    Author: Kyoung Tak Cho
    Created: Mon Jul 23 23:25:05 CDT 2018
    Updated: Tue Sep 08 21:41:58 CDT 2020
'''

import settings
import sqls
from utils.DBManager import Pgsql
from utils.Common import ListTool, StrTool, VersionManager
import random
import copy
import pprint
import Kmer as km

#
# GeneSequence class
#
class GeneSequence(object):
    def __init__(self, gnid=None, seq_type='', is_max_seq_len=None, conn=None, k=3):
        '''
        Constructor
        :param seq_str:
            sequence string
        :param seq_type:
            sequence types: string
                'p': Peptide sequence
                'd': DNA sequence
        '''
        if gnid is None:
            error_mesg = 'Error: gnid is null.'
            raise ValueError(error_mesg)
        self.gsid = None
        self.gnid = gnid
        self.seq_str = None
        self.seq_type = seq_type
        self.sub_id = None
        self.is_max_seq_len = is_max_seq_len
        self.k = k
        self.kmer_freq = [None] * settings.KMER_MAX_WINDOW_SIZE
        self._conn = conn

        #
        # Initialization
        #
        # set sequence information by gnid
        if self._set_seq_info() != 0:
            # k-mer freq build - default: 3-mer and 2-mer
            self._set_kmer_freq(k=k)
            self._set_kmer_freq(k=k-1)

    def len(self):
        '''
        Get sequence length
        :return:
            len(seq_str)
        '''
        return len(self.seq_str)

    def get_gsid(self):
        '''
        Get gsid
        :return:
            gsid: int
        '''
        return self.gsid

    def get_gnid(self):
        '''
        Get gnid if None, get from database
        :return:
            gnid: int
        '''
        if self.gsid is None:
            print('Error: no gsid in Seq.')
            return

        if self.gnid is None:
            res = Pgsql.Common.select_data(sqls.get_gnid_by_gsid, self.gsid)
            self.gnid = res[0]

        return self.gnid

    def get_seq_str(self):
        '''
        Get sequence string
        :return:
            seq_str: string
        '''
        return self.seq_str

    def _set_seq_info(self):
        '''
        Set gene sequence information from database
        :return:
            len(seq_str)
        '''

        # get sequence detail information from DB
        if self.is_max_seq_len:
            sql = sqls.get_seq_info_by_gnid_seq_len.format(MM='DESC')   # get max(seq_len)
        else:
            sql = sqls.get_seq_info_by_gnid_seq_len.format(MM='ASC')    # get min(seq_len)
        res = Pgsql.Common.select_data(sql=sql, pars=(self.gnid, self.seq_type), cur=self._conn.cursor())

        # set sequence detail information
        try:
            self.gsid = res[0][0]
            self.sub_id = res[0][1]
            self.seq_str = res[0][2]
        except IndexError as e:
            print(e)
            #print(len(res), res)
            print(sql % (self.gnid, self.seq_type))
            #print(self._conn)
            return 0

        if self.seq_str is None:
            return 0
        else:
            return (len(self.seq_str))

    def get_seq_type(self):
        '''
        Get sequence type
        :return:
            seq_type: string (char(1))
        '''
        return self.seq_type

    def set_seq_type(self, seq_type):
        '''
        Set sequence type
        :return:
            None
        '''
        self.seq_type = seq_type

    def lower(self):
        '''
        get lower case of sequence string
        :return:
            str(seq_str).lower()
        '''
        return str(self.seq_str).lower()

    def upper(self):
        '''
        get upper case of sequence string
        :return:
            str(seq_str).upper()
        '''
        return str(self.seq_str).upper()

    def get_kmer_freq(self, k=None):
        '''
        Get kmer frequency dictionary
            if kmer_freq is None:
                get kmer_freq from database
        :return:
            kmer_freq: dict()
        '''
        if k is None:
            error_mesg = 'k (window size) is empty.'
            raise ValueError(error_mesg)
        if self.seq_str is None:
            error_mesg = 'sequence string is empty.'
            raise ValueError(error_mesg)
        if self.kmer_freq[k-1] is None:
            # get k-mer frequency data from Database
            self._set_kmer_freq(k=k)

        return self.kmer_freq[k-1].get_kmer_freq()

    def _set_kmer_freq(self, k=None):
        if k is None:
            error_mesg = 'k (window size) is empty.'
            raise ValueError(error_mesg)

        kmer_freq = KmerFreq(seq=self, k=k, conn=self._conn)
        self.kmer_freq[k-1] = kmer_freq

#
# K-mer frequency class
#
class KmerFreq(object):
    def __init__(self, seq=None, k=None, conn=None):
        '''
        Constructor
            create KmerFreq object with sequence string and k (window size of k-mer)
        :param seq:
            sequence instance: GeneSequence
        :param k:
            window size k
        '''
        if seq is None:
            error_mesg = 'sequence instance is empty.'
            raise ValueError(error_mesg)
        if k is None:
            error_mesg = 'k (window size) is empty.'
            raise ValueError(error_mesg)

        self.k = k
        self.seq = seq
        self.kmer_freq = FreqDict(self.k)   # FreqDict
        self._conn = conn

        #
        # Initialization
        #

        # set Kmer Frequency
        if settings.RA_MODE:
            self._set_kmer_freq_ra()
        else:
            self._set_kmer_freq_dynm()

    def get_kmer_freq(self):
        if self.seq is None:
            error_mesg = 'sequence instance is empty.'
            raise ValueError(error_mesg)
        if self.kmer_freq is None:
            self._set_kmer_freq_dynm()
        return copy.copy(self.kmer_freq)

    def _set_kmer_freq(self, seq_type=None):
        if self.seq is None:
            error_mesg = 'sequence instance is empty.'
            raise ValueError(error_mesg)

        sql = sqls.get_kmer_freq_by_gsid
        if seq_type is None:        # default - peptide sequence
            pars = ('kmers', self.k, self.seq.get_gsid())
        elif seq_type == 'm1':       # for promoter sequence
            pars = ('kmers_promoter', self.k, self.seq.get_gsid())
        elif seq_type == 'd':       # for DNA sequence
            pars = ('kmers', self.k, self.seq.get_gsid())
        else:                       # default - peptide sequence
            pars = ('kmers', self.k, self.seq.get_gsid())

        if self._conn.closed:
            cur = None
        else:
            cur = self._conn.cursor()
        res = Pgsql.Common.select_data(sql=sql,
                                       pars=pars,
                                       cur=cur)
        for row in res:
            kmer = str(row[0]).upper()
            freq = row[1]
            self.kmer_freq.add_kmer_freq(kmer=kmer, freq=freq)

    def _set_kmer_freq_dynm(self):
        if self.seq is None:
            error_mesg = 'sequence instance is empty.'
            raise ValueError(error_mesg)

        kf = km.Kmer(seq=self.seq.get_seq_str(), window_size=self.k)
        #print('gsid: ', self.seq.get_gsid(), 'seq_len: ', self.seq.len(), self.seq.get_seq_str())
        for kmer, freq in kf.kmer_freq.items():
            self.kmer_freq.add_kmer_freq(kmer=kmer, freq=freq)
        #print('kmer_freq_len: ', len(self.kmer_freq))

    def print_kmer_freq(self):
        self.kmer_freq.print_freq()

    def set_by_gnid(self, gnid=None):
        pass


#
# Frequency Dictionary
#
class FreqDict(object):
    def __init__(self, k=None):
        self.k = k
        self.kmer_freq = dict()

    def __iter__(self):
        for kf in self.kmer_freq:
            yield kf

    def __len__(self):
        return len(self.kmer_freq)

    def __add__(self, other):
        '''
            Addition two KmerFreq object
        :param other:
        :return:
        '''
        added = FreqDict(k=self.k)
        # make a copy
        for kmer in self.kmer_freq:
            added.kmer_freq[kmer] = self.kmer_freq[kmer]
        # adding
        for kmer in other.kmer_freq:
            kmer = str(kmer).upper()
            added.kmer_freq[kmer] = self.kmer_freq.get(kmer, 0) + other.kmer_freq[kmer]
        return added

    def __sub__(self, other):
        sub = FreqDict(k=self.k)
        # make a copy
        for kmer in self.kmer_freq:
            sub.kmer_freq[kmer] = self.kmer_freq[kmer]
        # subtracting
        for kmer in other.kmer_freq:
            kmer = str(kmer).upper()
            cur_freq = sub.kmer_freq.get(kmer, None)
            if cur_freq is None:
                print('WARNING: key({}) missing in FreqDict.'.format(kmer))
            elif cur_freq > 0:
                new_freq = cur_freq - other.kmer_freq[kmer]
                if new_freq > 0:
                    # update frequency
                    sub.kmer_freq[kmer] = new_freq
                elif new_freq == 0:
                    del sub.kmer_freq[kmer]
                else:   # new_freq < 0
                    print('WARNING: new key({}) has negative value: {}'.format(kmer, new_freq))
                    del sub.kmer_freq[kmer]
            elif cur_freq == 0:
                print('WARNING: Zero-frequency value exists in in FreqDict. key:{}'.format(kmer))
                del sub.kmer_freq[kmer]
            else:   #cur_freq < 0:
                print('WARNING: current key({}) has negative value: {}'.format(kmer, cur_freq))
                del sub.kmer_freq[kmer]
        return sub

    def set_kmer_freq(self, kmer_freq):
        self.kmer_freq = kmer_freq

    def add_kmer_freq(self, kmer, freq):
        '''
            Add a kmer frequency into the dictionary
        :param kmer:
        :param freq:
        :return:
        '''
        kmer = str(kmer).upper()
        self.kmer_freq[kmer] = self.kmer_freq.get(kmer, 0) + freq

    def sub_kmer_freq(self, kmer, freq):
        new = copy.copy(self)
        kmer = str(kmer).upper()
        old_freq = new.kmer_freq.get(kmer, 0)
        if old_freq != 0:
            new.kmer_freq[kmer] = new.kmer_freq[kmer] - freq

    def get_kmer_freq_value(self, kmer=None):
        if kmer is None:
            error_mesg = 'kmer is empty.'
            raise ValueError(error_mesg)

        return self.kmer_freq.get(kmer, .000001)

    def total_freq(self):
        total = sum(self.kmer_freq.values())
        return total

    def total_freq_float(self):
        return float(self.total_freq())

    def print_freq(self, kmer_freq=None):
        if kmer_freq is None:
            kmer_freq = self.kmer_freq
        if not isinstance(kmer_freq, dict()):
            raise TypeError('given kmer_freq is not typeof dict()')

        for kmer, freq in kmer_freq.items():
            print(kmer, freq)

    def print(self, kmer_freq=None, sort_type=0, limit=0):
        if kmer_freq is None:
            kmer_freq = self.kmer_freq

        if not isinstance(kmer_freq, dict):
            raise TypeError('given kmer_freq is not type of dict().')

        if limit == 0 or limit > len(kmer_freq):
            limit = len(kmer_freq)

        if sort_type == 0:      # not sorted
            for kmer, freq in list(kmer_freq.items())[:limit]:
                print(kmer, freq)
        elif sort_type == 1:    # sorted by key - ASC
            for kmer in sorted(kmer_freq)[:limit]:
                print(kmer, kmer_freq[kmer])
        elif sort_type == 2:    # sorted by key - DESC
            for kmer in sorted(kmer_freq, reverse=True)[:limit]:
                print(kmer, kmer_freq[kmer])
        elif sort_type == 3:    # sorted by value - ASC
            for kmer, freq in sorted(kmer_freq.items(), key=lambda kv: (kv[1], kv[0]))[:limit]:
                print(kmer, freq)
        elif sort_type == 4:    # sorted by value - DESC
            for kmer, freq in sorted(kmer_freq.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)[:limit]:
                print(kmer, freq)

    def to_list(self, key_value=0, sort_type=0, limit=0):

        ret_list = list()

        if limit == 0 or limit > len(self.kmer_freq):
            limit = len(self.kmer_freq)

        if sort_type == 0:      # not sorted
            for kmer, freq in list(self.kmer_freq.items())[:limit]:
                if (key_value == 1):    # value
                    ret_list.append(freq)
                elif (key_value == 2):  # key
                    ret_list.append(kmer)
                else:                   # default - value
                    ret_list.append(freq)
        elif sort_type == 1:    # sorted by key - ASC
            for kmer in sorted(self.kmer_freq)[:limit]:
                if (key_value == 1):    # value
                    ret_list.append(self.kmer_freq[kmer])
                elif (key_value == 2):  # key
                    ret_list.append(kmer)
                else:                   # default - value
                    ret_list.append(self.kmer_freq[kmer])
        elif sort_type == 2:    # sorted by key - DESC
            for kmer in sorted(self.kmer_freq, reverse=True)[:limit]:
                if (key_value == 1):    # value
                    ret_list.append(self.kmer_freq[kmer])
                elif (key_value == 2):  # key
                    ret_list.append(kmer)
                else:                   # default - value
                    ret_list.append(self.kmer_freq[kmer])
        elif sort_type == 3:    # sorted by value - ASC
            for kmer, freq in sorted(self.kmer_freq.items(), key=lambda kv: (kv[1], kv[0]))[:limit]:
                if (key_value == 1):    # value
                    ret_list.append(freq)
                elif (key_value == 2):  # key
                    ret_list.append(kmer)
                else:                   # default - value
                    ret_list.append(freq)
        elif sort_type == 4:    # sorted by value - DESC
            for kmer, freq in sorted(self.kmer_freq.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)[:limit]:
                if (key_value == 1):    # value
                    ret_list.append(freq)
                elif (key_value == 2):  # key
                    ret_list.append(kmer)
                else:                   # default - value
                    ret_list.append(freq)

        return ret_list




#
# Validation data group
#   - for single class and it has a group per folder
#   - CrossValidation support multi classes
#
class ValidationDataset(object):
    def __init__(self, gnids=None, fold_size=None):
        if gnids is None:
            error_mesg = 'gene ids is empty.'
            raise ValueError(error_mesg)
        if fold_size is None:
            error_mesg = 'fold size is empty.'
            raise ValueError(error_mesg)

        self.gnids = gnids
        self.fold_size = fold_size
        self.group = None   # gnids * fold_size
        self.confusion_matrix = None    # ConfusionMatrix * fold_size

        #
        # Initialization
        #
        # random partition
        self._random_partition()

    def __len__(self):
        return len(self.gnids)

    def __iter__(self):
        for fold in self.group:
            yield fold

    def get_fold_size(self):
        return len(self.group)

    def _random_partition(self):
        random.shuffle(self.gnids)
        # initialize group
        self.group = []
        self.confusion_matrix = []
        for i in range(0, self.fold_size):
            self.group.append([])
            # initialize confusion matrix
            self.confusion_matrix.append(ConfusionMatrix())
        for idx, gnid in enumerate(self.gnids):
            self.group[idx % self.fold_size].append(gnid)

    def set_confusion_matrix(self, cm=None, fold_idx=None):
        if cm is None:
            error_mesg = 'cm (confusion matrix) is empty.'
            raise ValueError(error_mesg)
        if fold_idx is None:
            error_mesg = 'fold_idx is empty.'
            raise ValueError(error_mesg)
        self.confusion_matrix[fold_idx].set_confusion_matrix(cm=cm)

    def get_confusion_matrix(self, fold_idx=None):
        if fold_idx is None:
            error_mesg = 'fold_idx is empty.'
            raise ValueError(error_mesg)
        if fold_idx < 0 and fold_idx >= self.fold_size:
            error_mesg = 'invalid fold_idx.'
            raise ValueError(error_mesg)
        return self.confusion_matrix[fold_idx]

    def get_gnids(self):
        return self.gnids


#
# Gene class
#
class Gene(object):
    def __init__(self, gnid=None, gene_id=None, gene_name=None, gene_symbol=None,
                 pep_seq_min=None, pep_seq_max=None, dna_seq=None):
        '''
        Constructor
            Gene class
        :param gnid:
        :param gene_id:
        :param gene_name:
        :param gene_symbol:
        :param pep_seq:
            DNA sequence : GeneSequence
        :param dna_seq:
            DNA sequence : GeneSequence
        '''
        if gnid is None:
            print('Error: no gnid.')
            return

        self.gnid = gnid
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.gene_symbol = gene_symbol
        self.pep_seq_min = pep_seq_min      # peptide sequence id (gsid)
        self.pep_seq_max = pep_seq_max      # peptide sequence id (gsid)
        self.dna_seq = dna_seq              # dna sequence id (gsid)


#
# Gene List
#   handling genes using list()
#
class GeneList(object):
    def __init__(self, gene_list=None):
        if gene_list is None:
            self.gene_list = list()
        else:
            self.set_gene_list(gene_list)

    def merge_to(self, other=None, remove_duplicates=True):
        if other is not None:
            res = self.gene_list + other.get_gene_list()
            res.sort()
            if remove_duplicates:
                res = list(dict.fromkeys(res))    # remove duplicates using dict()
        else:
            print('Warning: other is None.')
            res = self.gene_list.copy()
        return GeneList(res)

    def add(self, gnid=None):
        if gnid is not None and gnid not in self.gene_list:
            self.gene_list.append(gnid)

    def remove(self, gnid=None):
        if gnid is not None and gnid in self.gene_list:
            self.gene_list.remove(gnid)

    def get_gene_list(self):
        self.gene_list.sort()
        return self.gene_list

    def set_gene_list(self, gene_list):
        src = gene_list.copy()
        src.sort()
        self.gene_list = src


#
# GeneDataSet
#   Gene dataset for classification including positive and negative classes
#
class GeneDataSet(object):
    def __init__(self, positive=None, negative=None):
        self.positive_class = GeneList()
        self.negative_class = GeneList()
        self.all_genes = GeneList()

        self.set_positive_class(positive)
        self.set_negative_class(negative)

    def set_positive_class(self, gene_list=None):
        if gene_list is not None and isinstance(gene_list, list):
            self.positive_class = GeneList(gene_list)

    def set_negative_class(self, gene_list=None):
        if gene_list is not None and isinstance(gene_list, list):
            self.negative_class = GeneList(gene_list)

    def merge(self):
        self.all_genes = self.positive_class.merge_to(self.negative_class)

    def get_positive_class(self):
        return self.positive_class

    def get_negative_class(self):
        return self.negative_class

    def get_all_genes(self):
        self.merge()
        return self.all_genes       # return GeneList() type

    def get_all_gnids(self):
        self.merge()
        return self.all_genes.get_gene_list()   # return list() type


#
# Gene Prediction Info class
#
class GenePrediction(object):
    def __init__(self, gnid=None, gsid=None, class_size=None, assigned_class=None, predicted_class=None):
        if gnid is None:
            error_mesg = 'gnid is empty.'
            raise ValueError(error_mesg)
        if class_size is None:
            error_mesg = 'class size is empty'
            raise ValueError(error_mesg)

        self.gnid = gnid
        self.gsid = gsid
        self.class_size = class_size
        self.assigned_class = assigned_class
        self.predicted_class = predicted_class
        self.log_p = list()

        #
        # Initialization
        #
        self._init_log_p()

    def _init_log_p(self):
        for class_num in range(0, self.class_size):
            self.log_p.append(0.)

    def get_assigned_predicted(self):
        return self.assigned_class, self.predicted_class

    def set_log_p(self, class_num=None, log_value=None):
        if class_num is None:
            error_mesg = 'class num is empty.'
            raise ValueError(error_mesg)
        self.log_p[class_num] = log_value

    def set_log_p_list(self, log_value_list):
        if self.class_size != len(log_value_list):
            error_mesg = 'log_value size is different from class_size.'
            raise ValueError(error_mesg)

        self.log_p = log_value_list

    def get_log_p(self, class_num=None):
        if class_num is None:
            error_mesg = 'class num is empty.'
            raise ValueError(error_mesg)
        return self.log_p[class_num]

    def set_predicted_class(self, predicted_class=None, log_p=None):
        if predicted_class is None:
            error_mesg = 'predicted class is empty.'
            raise ValueError(error_mesg)

        # set pridicted class
        self.predicted_class = predicted_class
        if log_p is not None:
            self.set_log_p_list(log_value_list=log_p)

    def get_predicted_class(self):
        return self.predicted_class

    def get_assigned_class(self):
        return self.assigned_class

    def to_list(self):
        res = list()
        res.append(self.gnid)
        res.append(self.gsid)
        res.append(self.assigned_class)
        res.append(self.predicted_class)
        for log_p in self.log_p:
            res.append(log_p)
        return res

    def get_header_list(self):
        header = list()

        header.append('gnid')
        header.append('Gene ID')
        header.append('gsid')
        header.append('Assigned Class')
        header.append('Predicted Class')

        # log(p) for each class
        for class_number in range(len(self.log_p)):
            header.append('log(p) - c{}'.format(class_number))

        return header

    def load_from_list(self, list_src):
        self.gnid = col[0]
        self.gsid = col[1]
        self.assigned_class = col[2]

#
# Confusion Matrix
#
class ConfusionMatrix(object):
    def __init__(self):
        self.tp = 0
        self.tn = 0
        self.fp = 0
        self.fn = 0

        # summary attributes
        self.total_ins = None
        self.total_c0 = None
        self.total_c1 = None
        self.precision = None
        self.recall = None
        self.fp_rate = None
        self.ppc_rate = None
        self.f_measure = None
        self.accuracy = None
        self.correct_ins = None
        self.incorrect_ins = None

    def __add__(self, other):
        new_cm = ConfusionMatrix()
        new_cm.set_tp(self.get_tp() + other.get_tp())
        new_cm.set_tn(self.get_tn() + other.get_tn())
        new_cm.set_fp(self.get_fp() + other.get_fp())
        new_cm.set_fn(self.get_fn() + other.get_fn())
        return new_cm

    def set_confusion_matrix(self, cm):
        if not isinstance(cm, ConfusionMatrix):
            error_mesg = 'cm is not an instance of ConfusionMatrix.'
            raise TypeError(error_mesg)
        # set values
        self.set_tp(cm.get_tp())
        self.set_tn(cm.get_tn())
        self.set_fp(cm.get_fp())
        self.set_fn(cm.get_fn())

    def get_tp(self):
        res = self.tp
        return res

    def get_tn(self):
        res = self.tn
        return res

    def get_fp(self):
        res = self.fp
        return res

    def get_fn(self):
        res = self.fn
        return res

    def set_tp(self, tp=0):
        self.tp = tp

    def set_tn(self, tn=0):
        self.tn = tn

    def set_fp(self, fp=0):
        self.fp = fp

    def set_fn(self, fn=0):
        self.fn = fn

    def add_tp(self):
        self.tp += 1

    def add_tn(self):
        self.tn += 1

    def add_fp(self):
        self.fp += 1

    def add_fn(self):
        self.fn += 1

    def get_total_c0(self):
        if self.total_c0 is None:
            self.total_c0 = self.tp + self.fn
        return self.total_c0

    def get_total_c1(self):
        if self.total_c1 is None:
            self.total_c1 = self.fp + self.tn
        return self.total_c1

    def get_total_ins(self):
        if self.total_ins is None:
            self.total_ins = self.get_total_c0() + self.get_total_c1()
        return self.total_ins

    def get_precision(self):
        if self.precision is None:
            self.precision = self._get_fraction(dividend=self.tp,
                                                divisor=self.tp + self.fp)
        return self.precision

    def get_recall(self):
        if self.recall is None:
            self.recall = self._get_fraction(dividend=self.tp,
                                             divisor=self.tp + self.fn)
        return self.recall

    def get_fp_rate(self):
        if self.fp_rate is None:
            self.fp_rate = self._get_fraction(dividend=self.fp,
                                              divisor=self.fp + self.tn)
        return self.fp_rate

    def get_ppc_rate(self):
        if self.ppc_rate is None:
            self.ppc_rate = self._get_fraction(dividend=self.tp + self.fp,
                                               divisor=self.get_total_ins())
        return self.ppc_rate

    def get_f_measure(self):
        if self.f_measure is None:
            tmp = self._get_fraction(self.get_precision() * self.get_recall(),
                                     self.get_precision() + self.get_recall())
            self.f_measure = 2 * tmp
        return self.f_measure

    def get_accuracy(self):
        if self.accuracy is None:
            self.accuracy = self._get_fraction(dividend=self.tp + self.tn,
                                               divisor=self.get_total_ins())
        return self.accuracy

    def get_correct_ins(self):
        if self.correct_ins is None:
            self.correct_ins = self.tp + self.tn
        return self.correct_ins

    def get_incorrect_ins(self):
        if self.incorrect_ins is None:
            self.incorrect_ins = self.fp + self.fn
        return self.incorrect_ins

    def _get_fraction(self, dividend, divisor):
        res = dividend / divisor if divisor != 0 else 0
        return res

    @staticmethod
    def get_header_list():
        header = list()

        # gene_set
        if settings.RA_MODE:
            header.append('Mapping ID')
        else:
            header.append('Gene set')
        # k-mer size
        header.append('K size')
        # Feature name
        header.append('Feature name')
        # Fold number (fold/overall)
        header.append('Fold')
        if settings.RA_MODE:
            header.append('Alphabet_size')

        # Total instances
        header.append('Total instances')
        # Confusion matrix
        header.append('True positive')
        header.append('False negative')
        header.append('False positive')
        header.append('True negative')
        # Correctly classified instances
        header.append('Correctly classified instances')
        # Incorrectly classified instances
        header.append('Incorrectly classified instances')
        # Accuracy
        header.append('Accuracy')
        # Precision
        header.append('Precision')
        # Recall
        header.append('Recall')
        # F-measure
        header.append('F-Measure')
        # FP rate
        header.append('FP rate')
        # PPC rate
        header.append('PPC rate')

        return header

    def to_list(self):
        summary = list()

        # Total instances
        summary.append(self.get_total_ins())
        # Confusion matrix
        summary.append(self.get_tp())
        summary.append(self.get_fn())
        summary.append(self.get_fp())
        summary.append(self.get_tn())
        # Correctly classified instances
        summary.append(self.get_correct_ins())
        # Incorrectly classified instances
        summary.append(self.get_incorrect_ins())
        # Accuracy
        summary.append(self.get_accuracy())
        # Precision
        summary.append(self.get_precision())
        # Recall
        summary.append(self.get_recall())
        # F-measure
        summary.append(self.get_f_measure())
        # FP rate
        summary.append(self.get_fp_rate())
        # PPC rate
        summary.append(self.get_ppc_rate())

        return summary




#
# Summary
#
class Summary(object):
    def __init__(self, cm):
        if not isinstance(cm, ConfusionMatrix):
            error_mesg = 'cm is not an instance of ConfusionMatrix.'
            raise TypeError(error_mesg)
        self.cm = cm

        # init attributes
        self.total_ins = None
        self.total_c0 = None
        self.total_c1 = None
        self.precision = None
        self.recall = None
        self.fp_rate = None
        self.ppc_rate = None
        self.f_measure = None
        self.accuracy = None
        self.correct_ins = None
        self.incorrect_ins = None

    def _init_attrs(self):
        self.total_ins = None
        self.total_c0 = None
        self.total_c1 = None
        self.precision = None
        self.recall = None
        self.fp_rate = None
        self.ppc_rate = None
        self.f_measure = None
        self.accuracy = None
        self.correct_ins = None
        self.incorrect_ins = None


    def get_total_c0(self):
        if self.total_c0 is None:
            self.total_c0 = self.cm.get_fp() + self.cm.get_tn()
        return self.total_c0

    def get_total_c1(self):
        if self.total_c1 is None:
            self.total_c1 = self.cm.get_tp() + self.cm.get_fn()
        return self.total_c1

    def get_total_ins(self):
        if self.total_ins is None:
            self.total_ins = self.get_total_c0() + self.get_total_c1()
        return self.total_ins

    def get_precision(self):
        if self.precision is None:
            self.precision = self._get_fraction(dividend=self.cm.get_tp(),
                                                divisor=self.cm.get_tp() + self.cm.get_fp())
        return self.precision

    def get_recall(self):
        if self.recall is None:
            self.recall = self._get_fraction(dividend=self.cm.get_tp(),
                                             divisor=self.cm.get_tp() + self.cm.get_fn())
        return self.recall

    def get_fp_rate(self):
        if self.fp_rate is None:
            self.fp_rate = self._get_fraction(dividend=self.cm.get_fp(),
                                              divisor=self.cm.get_fp() + self.cm.get_tn())
        return self.fp_rate

    def _get_fraction(self, dividend, divisor):
        res = dividend / divisor if divisor != 0 else 0
        return res

#
# Feature
#
class Feature(object):
    def __init__(self, name=None, class_size=2,
                 assigned_genes=None, corresp_tissue=None):

        if name is None:
            print('Error: feature name is empty.')
            return
        if assigned_genes is None:
            print('Error: assigned feature is empty.')
            return

        self.name = name
        self.class_size = class_size
        self.assigned_genes = assigned_genes
        self.corresp_tissue = corresp_tissue
        self.dataset = list()   # feature dataset for each class
        self.prediction_results = dict()
        self.cm_set = list()    # ConfusionMatrix * (fold_size + overall)

    def set_dataset(self, wd_all_gnids=None, fold_size=None):
        if wd_all_gnids is None:
            error_mesg = 'wd_all_gnids is empty.'
            raise ValueError(error_mesg)

        if fold_size is None:
            fold_size = 1
            print ('fold size is empty. set fold_size = 1 (default)')
        if fold_size <= 0:
            fold_size = 1
            print ('fold size is less than zero. set fold_size = 1 (default)')

        for class_num in range(0, self.class_size):
            self.dataset.append(class_num)
            self.dataset[class_num] = list()
            sub_dataset = list()
            if class_num == 0:
                sub_dataset = list(set(wd_all_gnids) - set(self.assigned_genes))
            elif class_num == 1:
                sub_dataset = self.assigned_genes
            self.dataset[class_num] = ValidationDataset(gnids=sub_dataset, fold_size=fold_size)

    def get_dataset_size(self):
        '''
        Get dataset size (list length) for each class
        :return:
            [length of each class]
        '''
        res = list()
        for each_class in self.dataset:
            print(len(each_class))
            res.append(len(each_class))
        return res

    def set_prediction_results(self, prediction_restuls=None):
        if prediction_restuls is None:
            error_mesg = 'prediction results is empty.'
            raise ValueError(error_mesg)
        self.prediction_results = prediction_restuls

    def get_confusion_matrix(self):
        return self.cm_set

    def set_confusion_matrix_set(self, cm_set=None):
        if cm_set is None:
            error_mesg = 'cm_set is empty.'
            raise ValueError(error_mesg)
        self.cm_set = cm_set

#
# GeneGroup class
#
class GeneGroup(object):
    def __init__(self, data, size):
        '''
        Constructor
            GeneGroup
        :param data:
        :param size:
        '''
        self.allgenes = data


#
# Feature Group
#
class SequenceGroup(object):
    def __init__(self, data, size):
        '''
        Constructor
        '''
        seq_group = data


#
# Feature infor
#
class FeatureInfo(object):
    def __init__(self, fsid=None, fs_name=None, gp_type=None,
                       description=None, class_size=None):
        self.fsid = fsid
        self.fs_name = fs_name
        self.gp_type = gp_type
        self.description = description
        self.class_size = class_size

    ###########################################################################
    # Set
    ###########################################################################
    def set_fsid(self, fsid=None):
        self.fsid = fsid

    def set_fs_name(self, fs_name=None):
        self.fs_name = fs_name

    def set_gp_type(self, gp_type):
        self.gp_type = gp_type

    def set_description(self, description):
        self.description = description

    def set_class_size(self, class_size):
        self.class_size = class_size

    ###########################################################################
    # Get
    ###########################################################################
    def get_fsid(self):
        return self.fsid

    def get_fs_name(self):
        return self.fs_name

    def get_gp_type(self):
        return self.gp_type

    def get_description(self):
        return self.description

    def get_class_size(self):
        return self.class_size


#
# Configuration
#
class Configuration(object):
    def __init__(self, fs_name=None, seq_type=None, target_features=None, gp_type=None, class_size=None):
        self.fs_info = FeatureInfo(fs_name=fs_name,
                                   gp_type=gp_type,
                                   class_size=class_size)
        self.seq_type = seq_type
        self.target_features = target_features
        self.fold_size = None
        self.kmer_size = None
        self.ignore_null = True
        self.gene_load_mode = None
        self.test_db = True
        self.test_mode = None
        self.debug_mode = None
        self.neg_class_mode = settings.NEG_CLASS_MODE_RND_S
        self.neg_class_set_size = 1
        self.assigned_genes_limit = [int(10) for x in range(1, 24)]  # default: 10 genes * 2 classes for 23 tissue
        self.seq_len = 0    # default: 0. if non-zero value then
        self.gene_dataset = list()      # GeneDataSet() * FEATURE_SIZE
        self.missing_gnids_in_promoter = None
        self.genes_info = None
        self.version = None

        #
        # Global data
        #
        # Cutoffs
        self.cutoffs = None

        #
        # Initialization
        #

        # create gene_dataset list() with FEATURE_SIZE
        self._init_gene_dataset()

    ###########################################################################
    # Init
    ###########################################################################
    def _init_gene_dataset(self, feature_size=settings.FEATURE_SIZE):
        for feature in range(1, feature_size + 1):
            self.gene_dataset.append(GeneDataSet())

    ###########################################################################
    # Set
    ###########################################################################
    def set_genes_info(self, genes_info):
        self.genes_info = genes_info

    def set_fs_info(self, fs_info=None):
        self.fs_info = fs_info

    def set_seq_type(self, seq_type=None):
        self.seq_type = seq_type

    def set_target_features(self, target_features):
        self.target_features = target_features

    def set_fold_size(self, fold_size):
        self.fold_size = fold_size

    def set_kmer_size(self, kmer_size):
        self.kmer_size = kmer_size

    def set_ignore_null(self, ignore_null):
        self.ignore_null = ignore_null

    def set_gene_load_mode(self, gene_load_mode):
        self.gene_load_mode = gene_load_mode

    def set_test_mode(self, test_mode):
        self.test_mode = test_mode

    def set_debug_mode(self, debug_mode):
        self.debug_mode = debug_mode

    def set_neg_class_mode(self, neg_class_mode):
        self.neg_class_mode = neg_class_mode

    def set_neg_class_set_size(self, neg_class_set_size):
        self.neg_class_set_size = neg_class_set_size

    # Set Cutoffs
    def set_cutoffs(self, cutoffs=None):
        self.cutoffs = cutoffs

    # set assigned genes limit
    def set_assigned_genes_limit(self, limit=None):
        if limit is None:
            self.assigned_genes_limit = [int(10) for x in range(1, 24)]  # default: 10 genes * 2 classes for 23 tissue
        else:
            self.assigned_genes_limit = limit

    # Set Gene DataSet:
    def set_gene_dataset(self, feature_id=None, positive_class=None, negative_class=None):
        if feature_id is None or feature_id <= 0 or feature_id > settings.FEATURE_SIZE:
            raise ValueError('featuer_id is empty.')

        if positive_class is not None:
            self.gene_dataset[feature_id - 1].set_positive_class(positive_class)
        if negative_class is not None:
            self.gene_dataset[feature_id - 1].set_negative_class(negative_class)

    # Set missing gnids in promoter sequence data
    def set_missing_gnids_in_promoter(self):
        sql = sqls.get_missing_gnids_in_promoter
        res = Pgsql.Common.select_data(sql=sql)
        self.missing_gnids_in_promoter = Pgsql.Common.to_list(res)

    ###################################
    # fs_info
    ###################################
    def set_fsid(self, fsid=None):
        self.fs_info.set_fsid(fsid)

    def set_fs_name(self, fs_name=None):
        self.fs_info.set_fs_name(fs_name)

    def set_gp_type(self, gp_type=None):
        self.fs_info.set_gp_type(gp_type)

    def set_class_size(self, class_size=None):
        self.fs_info.set_class_size(class_size)

    # Store Version information
    def set_version(self, version):
        self.version = VersionManager(version)

    def is_ignore_null(self):
        return self.ignore_null

    ###########################################################################
    # Get
    ###########################################################################
    def get_genes_info(self):
        return self.genes_info

    def get_fs_info(self):
        return self.fs_info

    def get_seq_type(self):
        return self.seq_type

    def get_target_features(self):
        return self.target_features

    def get_fold_size(self):
        return self.fold_size

    def get_kmer_size(self):
        return self.kmer_size

    def get_ignore_null(self):
        return self.ignore_null

    def get_gene_load_mode(self):
        return self.gene_load_mode

    def get_test_mode(self):
        return self.test_mode

    def get_debug_mode(self):
        return self.debug_mode

    def get_neg_class_mode(self):
        return self.neg_class_mode

    def get_neg_class_set_size(self):
        return self.neg_class_set_size

    # Get Cutoffs
    def get_cutoffs(self):
        return self.cutoffs

    # get assigned genes limit
    def get_assigned_genes_limit(self):
        return self.assigned_genes_limit

    # get gene dataset (GeneDataSet)
    def get_gene_dataset(self, feature_id=None):
        if feature_id is None or feature_id <= 0 or feature_id > settings.FEATURE_SIZE:
            raise ValueError('feature_id is empty.')
        return self.gene_dataset[feature_id - 1]

    # get gnids list per feature_id
    def get_gene_dataset_gnids_list(self, feature_id=None):
        return self.get_gene_dataset(feature_id=feature_id).get_all_gnids()

    # get gnids for positive class
    def get_gene_dataset_pos(self, feature_id=None):
        return self.get_gene_datset(feature_id=feature_id).get_positive_class().get_gene_list()

    # get gnids for negative class
    def get_gene_dataset_neg(self, feature_id=None):
        return self.get_gene_dataset(feature_id=feature_id).get_negative_class().get_gene_list()

    # get all gnids in dataset
    def get_gene_dataset_all_gnids_list(self):
        gene_list_all = GeneList()
        for gene_dataset in self.gene_dataset:
            gene_list = GeneList(gene_dataset.get_all_gnids())
            gene_list_all = gene_list.merge_to(gene_list_all)
        return gene_list_all.get_gene_list()

    # get missing gnids in promoter
    def get_missing_gnids_in_promoter(self):
        if self.missing_gnids_in_promoter is None:
            self.set_missing_gnids_in_promoter()
        return self.missing_gnids_in_promoter

    ###################################
    # fs_info
    ###################################
    def get_fsid(self):
        return self.fs_info.get_fsid()

    def get_fs_name(self):
        return self.fs_info.get_fs_name()

    def get_gp_type(self):
        return self.fs_info.get_gp_type()

    def get_class_size(self):
        return self.fs_info.get_class_size()

    # get version
    def get_version(self):
        return self.version


#
# Cutoff class - for single cutoff data
#
class Cutoff(object):
    def __init__(self, raw_data=None):
        if raw_data is not None:
            self.gp_type = raw_data[0]
            self.tissue_id = raw_data[1]
            self.percentile = raw_data[2]
            self.cutoff = raw_data[3]

    def gp_type(self):
        return self.gp_type()

    def tissue_id(self):
        return self.tissue_id()

    def percentile(self):
        return self.percentile()

    def cutoff(self):
        return self.cutoff()


#
# Cutoff data with dictionary
#
class Cutoffs(object):
    def __init__(self):
        self.data = dict()  # multi keys dictionary

    def __len__(self):
        return len(self.data)

    def __iter__(self):
        for cf in self.data:
            yield cf

    def _add(self, raw_data=None):
        if raw_data is not None:
            self.data.update({(raw_data[0], raw_data[1], raw_data[2]): raw_data[3]})

    def cutoff(self, cutoff_key):
        if cutoff_key not in self.data:
            # Get missing cutoff data from Database
            self.query_cutoff((int(float(cutoff_key[2]) * 100)))
            print(cutoff_key, 'added.')

        cutoff = self.data[cutoff_key]
        return cutoff

    def query_cutoffs(self, percentile_range=None):
        """
        :param gp_type:
        :param tissue_id:
        :param percentile_range:
            string type of target percentiles range
            ex) '70, 100, 5': 70%, 75%, ..., 95%
            ex) '95, 69, -5': 95%, 90%, ..., 70%
            ex) '30,  0, 5':
        :return:
        """
        if percentile_range is None:
            percentile_range = '95, 100, 5'  # default percentile is TOP 5%

        # set percentile range
        p_range = StrTool.str2range(percentile_range)
        percentiles = ", ".join('0.%02d' % x for x in p_range)
        sql = sqls.get_cutoffs_by_percentiles.format(PCTS=percentiles)
        rows = Pgsql.Common.select_data(sql)
        for row in rows:
            self._add(row)

    def query_cutoff(self, percentile=None):
        if percentile is None:
            raise ValueError('percentile is empty.')
        percentile_range = '{}, {}, 5'.format(percentile, percentile + 1)
        self.query_cutoffs(percentile_range)


# GE&PA comb data class
class CombConf(object):
    def __init__(self, fsid=None, fs_name=None,
                 gp_type=None,
                 hl1=None, hl2=None,
                 gp1=None, gp2=None,
                 pt1=None, pt2=None,
                 condition=None):
        self.fsid = fsid
        self.fs_name = fs_name
        self.gp_type = gp_type
        self.hl1 = hl1
        self.hl2 = hl2
        self.gp1 = gp1
        self.gp2 = gp2
        self.pt1 = pt1
        self.pt2 = pt2
        self.condition = condition


