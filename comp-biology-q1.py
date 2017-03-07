from __future__ import division
import re
import math

seqNeg = []
f = open('promoterNegative.fa', 'r')
for line in f:
    if line[0] != '>':
       seqNeg.append(line)

n_seqNeg = len(seqNeg)
_Ew = ((len(seqNeg[0]) - 1) - 6 + 1) * n_seqNeg

seq = []
f = open('promoterPositive.fa', 'r')
for line in f:
    if line[0] != '>':
       seq.append(line)

n_seq = len(seq)
print(n_seq)
print(seq)

nt = ['A', 'C', 'G', 'T', '[A|C]', '[A|G]', '[A|T]', '[C|G]', '[C|T]', '[G|T]', '[A|C|G|T]']
pnt = [0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1]
motif = ''
consensus_seq = ''
prev_consensus_seq = ''

max_zscore = 0
prev_max_zscore = 0
max_Nw = prev_max_Nw = 0
max_Ew = prev_max_Ew = 0

print("motif \t zscore \t Nw \t Ew")

for n1 in range(0, 11):
    for n2 in range(0, 11):
        for n3 in range(0, 11):
            for n4 in range(0, 11):
                for n5 in range(0, 11):
                    _motif = nt[n1] + nt[n2] + nt[n3] + nt[n4] + nt[n5]
                    for n6 in range(0, 11):
                        motif = _motif
                        motif += nt[n6]
                        Nw = 0
                        Ew = 0
                        zscore = 0
                        for i in range(0, n_seq):
                            matches = re.findall(motif, seq[i], 1)
                            if matches:
                                Nw += 1
                        Ew = _Ew * pnt[n1] * pnt[n2] * pnt[n3] * pnt [n4] * pnt[n5] * pnt[n6]
                        zscore = (Nw - Ew)/math.sqrt(Ew)
                        if zscore >= max_zscore:
                            prev_max_zscore = max_zscore
                            prev_consensus_seq = consensus_seq
                            prev_max_Ew = max_Ew
                            prev_max_Nw = max_Nw
                            max_zscore = zscore
                            consensus_seq = motif
                            max_Ew = Ew
                            max_Nw = Nw

print("{} \t {} \t {} \t {}".format(consensus_seq, max_zscore, max_Nw, max_Ew))
print("{} \t {} \t {} \t {}".format(prev_consensus_seq, prev_max_zscore, prev_max_Nw, prev_max_Ew))
