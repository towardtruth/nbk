'''
   Kmer class
      Author: Kyoung Tak Cho
      Created: Thu Jul 19 00:39:38 CDT 2018
      Updated: Sun Aug  5 23:58:33 CDT 2018
'''


class Kmer(object):
    """ Kmer class
        Descriptions:
            Kmer class is designed for kmer objects for a sequence.
        Attributes:
            seq
            kmer
        Methods:
                compute kmer from a sequence
    """

    def __init__(self, seq=None, window_size=3, kmer_freq=None):
        if kmer_freq is None:
            if seq is None:
                raise ValueError('sequence is empty.')
            self.kmer_freq = Kmer.kmer_freq_acc(seq=seq, window_size=window_size)
        else:
            if seq is None:
                self.kmer_freq = kmer_freq
            else:
                self.kmer_freq = Kmer.kmer_freq_acc(seq=seq, window_size=window_size, kmer_freq=kmer_freq)

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

    @staticmethod
    def kmer_freq_acc(seq, window_size, freq=None):
        seq = str(seq).upper()

        # if freq is not dictionary type, makes it as a dict() type
        if freq is None or not isinstance(freq, dict):
            freq = dict()

        # compute all kmers frequency
        for i in range(0, len(seq) - (window_size - 1)):
            s = seq[i:i+window_size]
            #if s.find('N') < 0:
            #    freq[s] = freq.get(s, 0) + 1
            freq[s] = freq.get(s, 0) + 1
        return freq

    @staticmethod
    def sort(kmer_freq=None, sort_type=0, limit=0):
        if kmer_freq is None:
            raise ValueError('kmer_freq is empty.')

        if limit == 0 or limit > len(kmer_freq):
            limit = len(kmer_freq)

        new_kmer_freq = dict()
        if sort_type == 0:      # not sorted
            for kmer, freq in list(kmer_freq.items())[:limit]:
                new_kmer_freq[kmer] = freq
                #print(kmer, freq)
        elif sort_type == 1:    # sorted by key - ASC
            for kmer in sorted(kmer_freq)[:limit]:
                new_kmer_freq[kmer] = kmer_freq[kmer]
                #print(kmer, kmer_freq[kmer])
        elif sort_type == 2:    # sorted by key - DESC
            for kmer in sorted(kmer_freq, reverse=True)[:limit]:
                new_kmer_freq[kmer] = kmer_freq[kmer]
                #print(kmer, kmer_freq[kmer])
        elif sort_type == 3:    # sorted by value - ASC
            for kmer, freq in sorted(kmer_freq.items(), key=lambda kv: (kv[1], kv[0]))[:limit]:
                new_kmer_freq[kmer] = freq
                #print(kmer, freq)
        elif sort_type == 4:    # sorted by value - DESC
            for kmer, freq in sorted(kmer_freq.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)[:limit]:
                new_kmer_freq[kmer] = freq
                #print(kmer, freq)
        return Kmer(new_kmer_freq)