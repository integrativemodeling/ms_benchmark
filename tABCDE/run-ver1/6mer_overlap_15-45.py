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
	      for i in range(0,5):
		     for j in range(0,5):
                            DA = sqrt(pow(b[0][i]-b[0][j], 2) + pow(b[1][i]-b[1][j], 2) + pow(b[2][i]-b[2][j], 2))
              #DA3A1 = sqrt(pow(b[0][2]-b[0][0], 2) + pow(b[1][2]-b[1][0], 2) + pow(b[2][2]-b[2][0], 2))
              #DA3A2 = sqrt(pow(b[0][2]-b[0][1], 2) + pow(b[1][2]-b[1][1], 2) + pow(b[2][2]-b[2][1], 2))
              #DA3B1 = sqrt(pow(b[0][2]-b[0][3], 2) + pow(b[1][2]-b[1][3], 2) + pow(b[2][2]-b[2][3], 2))
              #DA3B2 = sqrt(pow(b[0][2]-b[0][4], 2) + pow(b[1][2]-b[1][4], 2) + pow(b[2][2]-b[2][4], 2))
              #DA3B3 = sqrt(pow(b[0][2]-b[0][5], 2) + pow(b[1][2]-b[1][5], 2) + pow(b[2][2]-b[2][5], 2))  
              return DA
	      #kl=getDistancesA3toALL(lines)
              #print kl
              

 
       def getOverlapsA3toALL(lines):
              for i in range(0,5):
                     for j in range(0,5):
                            OvA = (b[3][i]+b[3][j]-c)/(b[3][i]+b[3][j])
              #OvA3A1 = (b[3][2]+b[3][0]-c[0])/(b[3][2]+b[3][0])
              #OvA3A2 = (b[3][2]+b[3][1]-c[1])/(b[3][2]+b[3][1])
              #OvA3B1 = (b[3][2]+b[3][3]-c[2])/(b[3][2]+b[3][3])
              #OvA3B2 = (b[3][2]+b[3][4]-c[3])/(b[3][2]+b[3][4])
              #OvA3B3 = (b[3][2]+b[3][5]-c[4])/(b[3][2]+b[3][5])

              return OvA


       c = getDistancesA3toALL(ReadFile(InputFileName))
       d= getOverlapsA3toALL(ReadFile(InputFileName))
       print d            
       
       
       if (d[0]<0.45) and (d[1]<0.45) and (d[2]<0.45) and (d[3]<0.45) and (d[4]<0.45):
	       if (d[0]>0.15) or (d[1]>0.15) or (d[2]>0.15) or (d[3]>0.15) or (d[4]>0.15):
                       Fname = 'ov15-45_test.'+str(i)+'.mfj'
                       f1 = open(Fname, 'w')
                       f1.write(Fname+'\n')
                       f1.write('1'+'\n')
                       f1.write('6'+'\n')
                       f1.write('ang'+'\n')
                       f1.write('none'+'\n')
                       f1.write('1.0000'+'\n')
                       for i in ReadFile(InputFileName)[6:]:
                               l = re.findall ("(\S+)", i) 
                               x = float(l[0])
                               y = float(l[1])
                               z = float(l[2])
                               r = float(l[3])
                               f1.write("%s%8.3f%s%8.3f%s%8.3f%7.1f%s" % ('    ',x,'      ',y,'      ',z,r,'      .00000')+'\n')
                       f1.close()
                                                                                                     
       

        
        
        
