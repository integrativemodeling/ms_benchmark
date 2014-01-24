#!/usr/bin/python

# code to implement the xlinking scoring for the proteasome module 1 (Rpn3/Rpn7/SEM1)

import re
import os
from math import sqrt, pow

#NAstrs = raw_input("enter the first number of structures: ")
#NBstrs = raw_input("enter the last number of structures: ")
#Ns = raw_input("enetr number of subunits: ")

#NA = int(NAstrs)
#NB = int(NBstrs)
#Nsub = int(Ns)
NA = 0
NB = 10000
Nsub = 7

list0 = []
list1 = []
list2 = []
list3 = []
list4 = []
list5 = []
list6 = []
list7 = []
listNames = []
# Read input file
for i in range(NA, NB):
    InputFileName = "output/mfj/configuration." + str(i) + ".mfj"
    if not os.path.exists(InputFileName):
        continue

    def ReadFile(InputFile):
        InputFile = open(InputFileName, 'r')  # Open file to read
        # string of lines from pdb file
        InputFileLines = InputFile.read().splitlines()
        InputFile.close()
        return InputFileLines

    def coordinates(lines):
        outX = []
        outY = []
        outZ = []
        outR = []
        for i in lines[6:]:
            l = re.findall("(\S+)", i)
            X0 = float(l[0])
            Y0 = float(l[1])
            Z0 = float(l[2])
            R0 = float(l[3])
            outX.append(X0)
            outY.append(Y0)
            outZ.append(Z0)
            outR.append(R0)
        return outX, outY, outZ, outR
    b = coordinates(ReadFile(InputFileName))
    # we dont care for the internal connections as those are fixed through distance restraints
    # As are Rpn3, Bs Rpn7 and Cs Sem1
    # distances between the cross linked residues nas6 (210) and rpt1 (220) - in models from IMP where structures
    # are coarse grained to residue level are identified as 208 and 493 (in
    # python: 207 and 492 tespectuively)

    def getDiffs(linkA, linkB, RA, RB):
        Diff1 = (
            sqrt(pow(b[0][linkA] - b[0][linkB],
                     2) + pow(b[1][linkA] - b[1][linkB],
                              2) + pow(b[2][linkA] - b[2][linkB],
                                       2))) - (b[3][RA] + b[3][RB]) * 1.25
        return Diff1

    # X-link interactions observed in experiments
    # module 1: -A2B3 -A3B3 -B3C2 -A1C2
    # SIDs:
    # 11-5 / E1-A3 >100 then penalty score=0.0
    # 11-8 / E1-C1 (>100) and E2-C2 (<50, >20)
    # 6-8 / B1 - C2 (>100)
    # 5 - 8 / (A3-C2) (20<x<50)
    # 5- 9 / A2-D1 (>100)
    # 9-8 / D1-C1 (>100)

    def getScoreInteractions(cLink):
        if cLink < 0:
            score1 = 0.0
        else:
            score1 = 2.0
        Score_Inter = score1
        return Score_Inter

    def getScoreNonInteractions(cLink):
        if cLink > 0.0:
            score2 = 0.0
        else:
            score2 = 1.0
        Score_nonInter = score2
        return Score_nonInter

    listNames.append(InputFileName)
    outDiff = []
    for i in range(0, Nsub):
        for j in range(0, Nsub):
            if i < j:
                Diff = getDiffs(i, j, i, j)
                outDiff.append(Diff)
    #print outDiff

    list1 = [0, 1, 5, 7, 12, 20]
    ScoreI = []
    for i in list1:
        SA = getScoreInteractions(outDiff[i])
        ScoreI.append(SA)
    SCI = sum(ScoreI)
    print SCI

    list2 = [2, 3, 4, 6, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19]
    ScoreNI = []
    for i in list2:
        SB = getScoreNonInteractions(outDiff[i])
        ScoreNI.append(SB)
    SCNI = sum(ScoreNI)
    print SCNI

    STot = SCI + SCNI
    list0.append(str(STot))

SummaryFileName = 'output/score-connect.txt'

SummaryFile = open(SummaryFileName, 'w')
SummaryFile.write(" " .join(listNames) + '\n')
SummaryFile.write(" ".join(list0) + '\n')
