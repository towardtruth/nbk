"""
    kmer_test.py
    author: Kyoung Tak Cho
    created: 03/09/2019
"""

import unittest
#from nbk.Kmer import Kmer
#from .. Kmer import Kmer
from .. import Kmer
from nbk.models import FreqDict


class UTestKmerFreq(unittest.TestCase):
    def test_kmer_freq(self):
        test_seq = 'ABCDEFG'
        test_window_size = 3

        result_freq = Kmer.kmer_freq_acc(test_seq, test_window_size)
        fd = FreqDict(3)
        fd.set_kmer_freq(result_freq)

        fd += fd

        print(result_freq)
        print(len(fd))


if __name__ == '__main__':
    unittest.main()

