from __future__ import division
import numpy as np
import re

class CrossValidation(object):
    '''
    Iterator that returns 1/k of the data as validation data and
    the rest as training data, for every of the k pieces.
    '''
    def __init__(self, examples, k=1):

        self.examples = examples
        self.k = k
        self.i = 0

    def __iter__(self):
        return self

    def next(self):
        s, e = self.i * self.k, (self.i + 1) * self.k
        if s >= len(self.examples):
            raise StopIteration
        self.i += 1

        x1 = self.examples[:s]
        x2 = self.examples[e:]
        train_data = np.concatenate((x1,x2))

        test_data = self.examples[s:e]

        return train_data, test_data

class Training(object):

     def __init__(self, train_seq, length, nPos):
        self.train_seq = train_seq
        self.length = length
        self.nPos = nPos
        self.nt = ['A', 'C', 'G', 'T', '[A|C]', '[A|G]', '[A|T]', '[C|G]', '[C|T]', '[G|T]', '[A|C|G|T]']
        self.consensus_seq = []
        self.seq = []
        self._seq = np.empty([self.length, self.nPos], dtype = str)

        for i in range(0, self.length):
            for j in range(0, self.nPos):
                self._seq[i, j] = self.train_seq[j][i]
            self.seq.append(''.join(self._seq[i,:]))

     def trainCS(self):
        for i in range(0, self.length):
            for j in range(0, len(self.nt)):
                matches = re.findall(self.nt[j], self.seq[i], 0)
                if float(len(matches)/self.nPos) >= 0.9:
                    self.consensus_seq.append(self.nt[j])
                    break

        return self.consensus_seq

     def trainPWM(self):
        pwm = np.empty([4, self.length], dtype = float)
        for j in range(0, self.length):
            for i in range(0, 4):
                matches = re.findall(self.nt[i], self.seq[j], 0)
                nmatches = len(matches)
                pwm[i,j] = float(nmatches / self.nPos)

        return pwm

     def trainPWMPseudoCnts(self):
        pwm = np.empty([4, self.length], dtype = float)
        for j in range(0, self.length):
            for i in range(0, 4):
                matches = re.findall(self.nt[i], self.seq[j], 0)
                nmatches = len(matches)
                pwm[i,j] = float((nmatches + 1) / (self.nPos + 4))

        return pwm
