from __future__ import division
import re
from random import uniform
import numpy as np
import utils
import matplotlib.pyplot as plt

def frange(start, end, step):
    tmp = start
    while(tmp < end):
        yield tmp
        tmp += step

def calculateLikelihood(j, length, sequences, pwm):
    likelihood = 1
    for i in range(0, length):
        if sequences[j][i] == 'A':
            likelihood *= pwm[0, i]
        elif sequences[j][i] == 'C':
            likelihood *= pwm[1, i]
        elif sequences[j][i] == 'G':
            likelihood *= pwm[2, i]
        elif sequences[j][i] == 'T':
            likelihood *= pwm[3, i]
    return likelihood

def loadDataset(filename):
    f = open(filename, 'r')
    sequences = []
    for line in f:
        sequences.append(line)
    return sequences

def drawFigure(title1, title2, label1, label2, data1, data2):
    plt.plot(data1, data2, linewidth=3.0)
    
    plt.xlabel(label1)
    plt.ylabel(label2)

    plt.title(title1 + ' on ' + title2)

    plt.savefig(title1 + 'on' + title2 + '.png')

    plt.show()
    
#-------------------------------------------------------------------------
dataset = 'dataset2'
negative = loadDataset('negative2.txt')
positive = loadDataset('positive2.txt')

stepSize = 0.0001

TruePos = TrueNeg = 0

nPos = len(positive) 
nNeg = len(negative) 
length = len(positive[0]) - 1
print("#of positive: {}".format(nPos))
print("#of negative: {}".format(nNeg))
print("length: {}".format(length))
#-------------------------------------------------------------------------
#LOOCV Consensus sequences
cv = utils.CrossValidation(positive)

for i in range(0, nPos):
    train_seq, test_seq = cv.next()
    t1 = utils.Training(train_seq, length, nPos - 1)
    consensus_seq = t1.trainCS()
    nmatch = 0
    for j in range(0, length):
        matches = re.findall(consensus_seq[j], test_seq[0][j], 0)
        if matches:
            nmatch += 1
    if nmatch == length:
        TruePos += 1
#Train all
t1 = utils.Training(positive, length, nPos)
consensus_seq = t1.trainCS()
print("consensus seq: {}".format(consensus_seq))
#True Negative
for i in range(0, nNeg):
    nmatch = 0
    for j in range(0, length):
        matches = re.findall(consensus_seq[j], negative[i][j], 0)
        if matches:
            nmatch += 1
    if nmatch != length:
        TrueNeg += 1
#Sensitivity and specificity
print("Consensus Sequence:")
print("Sensitivity = TP/n = {}/{} = {}".format(TruePos, nPos, float(TruePos/nPos)))
print("Specificity = TN/m = {}/{} = {}".format(TrueNeg, nNeg, float(TrueNeg/nNeg)))
#-------------------------------------------------------------------------
#LOOCV PWM
sensitivity = []
specificity = []
FPFN = []
t = []

threshold = frange(0, 1, stepSize)

for thrsh in threshold:
    FalsePos = FalseNeg = 0
    TruePos = TrueNeg = 0
    cv = utils.CrossValidation(positive)
    for j in range(0, nPos):
        train_seq, test_seq = cv.next()
        t1 = utils.Training(train_seq, length, nPos - 1)
        pwm = t1.trainPWM()
        likelihood = calculateLikelihood(0, length, test_seq, pwm)
        if likelihood > thrsh:
            TruePos += 1
        else:
            FalseNeg += 1

    #Train all
    t1 = utils.Training(positive, length, nPos)
    pwm = t1.trainPWM()
    #True Negative
    for j in range(0, nNeg):
        likelihood = calculateLikelihood(j, length, negative, pwm)
        if likelihood <= thrsh:
            TrueNeg += 1
        else:
            FalsePos += 1

    #Sensitivity and specificity
    sensitivity.append(float(TruePos/nPos))
    specificity.append(float(TrueNeg/nNeg))
    FPFN.append(FalsePos + FalseNeg)
    t.append(thrsh)

print("PWM")
print(pwm)
idx = np.argmin(FPFN)
print('FalsePos + FalseNeg: {}'.format(FPFN[idx]))
print('Threshold: {}'.format(t[idx]))
print('Sensitivity: {}'.format(sensitivity[idx]))
print('Specificity: {}'.format(specificity[idx]))

drawFigure('Threshold', dataset, 'Threshold', 'False-positives + False-negatives', t, FPFN)

drawFigure('PWM', dataset, 'Specificity', 'Sensitivity', specificity, sensitivity)
#-------------------------------------------------------------------------
#LOOCV PWM
sensitivity = []
specificity = []

threshold = frange(0, 1, stepSize)

for thrsh in threshold:
    TruePos = TrueNeg = 0
    cv = utils.CrossValidation(positive)
    for j in range(0, nPos):
        train_seq, test_seq = cv.next()
        t2 = utils.Training(train_seq, length, len(train_seq) - 1)
        pwm = t2.trainPWMPseudoCnts()
        likelihood = calculateLikelihood(0, length, test_seq, pwm)
        if likelihood > thrsh:
            TruePos += 1

    #Train all
    t1 = utils.Training(positive, length, nPos)
    pwm = t1.trainPWMPseudoCnts()
    #True Negative
    for j in range(0, nNeg):
        likelihood = calculateLikelihood(j, length, negative, pwm)
        if likelihood <= thrsh:
            TrueNeg += 1

    #Sensitivity and specificity
    sensitivity.append(float(TruePos/nPos))
    specificity.append(float(TrueNeg/nNeg))

print("Pseudocounts PWM")
print(pwm)

drawFigure('Pseudocounts PWM', dataset, 'Specificity', 'Sensitivity', specificity, sensitivity)
