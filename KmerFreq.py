'''
    Kmer Frequency class
    Author: Kyoung Tak Cho
    Created: June 4, 2018
'''


class KmerFreq(dict):

    #def __init__(self):

    def get_kmer_freq(seq, window_size, freqAll):
        seq = str(seq).upper()

        for i in range(0, len(seq) - (window_size - 1)):
            s = seq[i:i + window_size]
            freqAll[s] = freqAll.get(s, 0) + 1

        return freqAll

    def get_kmer_freq(seq, window_size):
        freq = dict()
        return get_ngram_freq(seq, window_size, freq)


