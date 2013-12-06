#import wx
import re
#import os, numarray
import operator
from operator import itemgetter
#from numpy  import *
#from Numeric import *
from math import *
#from Numeric import * # imports numerical python
import csv


Nstart = raw_input("What is start number of models?:  ")
Nend = raw_input("What is end+1 number of models?:  ")
Na = int(Nstart)
Nb = int(Nend)


list0 = []
listNames=[]
##Read input file
for i in range(Na, Nb):
       InputFileName="mod-abcde-ver1."+str(i)+".mfj"
       def ReadFile(InputFile):
              InputFile = open(InputFileName,'r') # Open file to read
              InputFileLines = InputFile.read().splitlines() # string of lines from pdb file
              InputFile.close()
              return InputFileLines


       def coordinates(lines):
              outX=[]
              outY=[]
              outZ=[]
              outR=[]
              for i in lines[6:]:
                     l = re.findall ("(\S+)", i) 
                     X0 = float(l[0])
                     Y0 = float(l[1])
                     Z0 = float(l[2])
                     R0 = float(l[3])
                     outX.append(X0)
                     outY.append(Y0)
                     outZ.append(Z0)
                     outR.append(R0)
                     #print Scr2
              return outX, outY, outZ, outR
       b= coordinates(ReadFile(InputFileName))
       #print b[0][0]
       coordinates(ReadFile(InputFileName))
       def getDistancesA3toALL(lines):

              DAB = sqrt(pow(b[0][0]-b[0][1], 2) + pow(b[1][0]-b[1][1], 2) + pow(b[2][0]-b[2][1], 2))
              DAC = sqrt(pow(b[0][0]-b[0][2], 2) + pow(b[1][0]-b[1][3], 2) + pow(b[2][0]-b[2][2], 2))
              DAD = sqrt(pow(b[0][0]-b[0][3], 2) + pow(b[1][0]-b[1][3], 2) + pow(b[2][0]-b[2][3], 2))
              DAE = sqrt(pow(b[0][0]-b[0][4], 2) + pow(b[1][0]-b[1][4], 2) + pow(b[2][0]-b[2][4], 2))
              DBC = sqrt(pow(b[0][1]-b[0][2], 2) + pow(b[1][1]-b[1][2], 2) + pow(b[2][1]-b[2][2], 2))  
              DBD = sqrt(pow(b[0][1]-b[0][3], 2) + pow(b[1][1]-b[1][3], 2) + pow(b[2][1]-b[2][3], 2))
              DBE = sqrt(pow(b[0][1]-b[0][4], 2) + pow(b[1][1]-b[1][4], 2) + pow(b[2][1]-b[2][4], 2))
              DCD = sqrt(pow(b[0][2]-b[0][3], 2) + pow(b[1][2]-b[1][3], 2) + pow(b[2][2]-b[2][3], 2))
              DCE = sqrt(pow(b[0][2]-b[0][4], 2) + pow(b[1][2]-b[1][4], 2) + pow(b[2][2]-b[2][4], 2))
              DDE = sqrt(pow(b[0][3]-b[0][4], 2) + pow(b[1][3]-b[1][4], 2) + pow(b[2][3]-b[2][4], 2))
              return DAB, DAC, DAD, DAE, DBC, DBD, DBE, DCD, DCE, DDE


       def getOverlapsA3toALL(lines):
              OvAB = (b[3][0]+b[3][1]-c[0])/(b[3][0]+b[3][1])
              OvAC = (b[3][0]+b[3][2]-c[1])/(b[3][0]+b[3][2])
              OvAD = (b[3][0]+b[3][3]-c[2])/(b[3][0]+b[3][3])
              OvAE = (b[3][0]+b[3][4]-c[3])/(b[3][0]+b[3][4])
              OvBC = (b[3][1]+b[3][2]-c[4])/(b[3][1]+b[3][2])
              OvBD = (b[3][1]+b[3][3]-c[5])/(b[3][1]+b[3][3])
              OvBE = (b[3][1]+b[3][4]-c[6])/(b[3][1]+b[3][4])
              OvCD = (b[3][2]+b[3][3]-c[7])/(b[3][2]+b[3][3])
              OvCE = (b[3][2]+b[3][4]-c[8])/(b[3][2]+b[3][4])
              OvDE = (b[3][3]+b[3][4]-c[9])/(b[3][3]+b[3][4])
              return OvAB, OvAC, OvAD, OvAE, OvBC, OvBD, OvBE, OvCD, OvCE, OvDE


       c = getDistancesA3toALL(ReadFile(InputFileName))
       d= getOverlapsA3toALL(ReadFile(InputFileName))
                  

       def getScore(lines):
              if (d[0]<0.45) and (d[1]<0.45) and (d[2]<0.45) and (d[3]<0.45) and (d[4]<0.45) and (d[5]<0.45) and (d[6]<0.45) and (d[7]<0.45) and (d[8]<0.45) and (d[9]<0.45):
                     score1 = 0.0
              else:
                     score1 = 1.0
              ScoreTotal = score1
              return ScoreTotal
       
       listNames.append(InputFileName)
       STot = getScore(ReadFile(InputFileName))
       eST=str(STot)





       list0.append(eST)

       SummaryFileName = 'out_overlap-upd.csv'
       SummaryFile = open(SummaryFileName, 'w')
       SummaryFile.write(str(listNames) +'\n')
#       SummaryFile.write(str(list0) +'\n')
#       SummaryFile.write(str(list1) +'\n')
 #      SummaryFile.write(str(list2) +'\n')
 #      SummaryFile.write(str(list3) +'\n')
 #      SummaryFile.write(str(list4) +'\n')
 #      SummaryFile.write(str(list5) +'\n')
 #      SummaryFile.write(str(list6) +'\n')
 #      SummaryFile.write(str(list7) +'\n')
       SummaryFile.write(str(list0) +'\n')


       #SummaryFile.write(str(listOv4) +'\n')
       #SummaryFile.write(str(listOv5) +'\n')
       #SummaryFile.write(str(listOv6) +'\n')


SummaryFile.flush()
SummaryFile.close()







#      if (d[0]>0.15) or (d[1]>0.15) or (d[2]>0.15) or (d[3]>0.15) or (d[4]>0.15) and (d[5]>0.15) or (d[6]>0.15) or (d[7]>0.15) or (d[8]>0.15) or (d[9]>0.15):
#                       Fname = 'ov-max45_modABCDE'+str(i)+'.mfj'
#                       f1 = open(Fname, 'w')
#                       f1.write(Fname+'\n')
#                       f1.write('1'+'\n')
#                       f1.write('6'+'\n')
#                       f1.write('ang'+'\n')
#                       f1.write('none'+'\n')
#                       f1.write('1.0000'+'\n')
#                       for i in ReadFile(InputFileName)[6:]:
#                               l = re.findall ("(\S+)", i) 
#                               x = float(l[0])
#                               y = float(l[1])
#                               z = float(l[2])
#                               r = float(l[3])
#                               f1.write("%s%8.3f%s%8.3f%s%8.3f%7.1f%s" % ('    ',x,'      ',y,'      ',z,r,'      .00000')+'\n')
#                       f1.close()
#                                                                                                     
       

        
        
        
