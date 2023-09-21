
from  genomeScan_lib import *


########################################################################
########################################################################
########################################################################
########################################################################


###Global Vars

##Stop condon are three TAG TAA TGA and some alternatives (to check if exist for human)
startCodonList=["AAG","ACG","AGG","ATA","ATC","ATG","ATT","CTG","GTG","TTG"]
startCodonList=["ATG"]

stopCodonList=["TAG","TAA","TGA"]
stopCodonList=["TGA"]

chromList=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
chromList=["22"]

startCodonPosList=[]
stopCodonPosList=[]

dictFileName="data/dnaDict.pkl"

#######################MAIN################################################
import sys

argv = sys.argv[1:]

if len (argv)==0:
    standAloneFlag=True
else:
    standAloneFlag=False
    aChrom = argv[0]
    aFrame = argv[1]
print ("run for chrom ",aChrom )
# timeStart()
dnaDict=readDict(dictFileName)
# timeEnd()

timeStart()

if standAloneFlag:
    aFrame=1
    for aChrom in chromList:
        startCodonPosList,stopCodonPosList = runFunction(aChrom,dnaDict,aFrame)
        print (aChrom,len(startCodonPosList),len(stopCodonPosList))
        findSorfsCoords(aChrom,startCodonPosList,stopCodonPosList)
    timeEnd()

if (standAloneFlag==False):
    startCodonPosList,stopCodonPosList = runFunction(aChrom,dnaDict[aChrom],aFrame)
    print (aChrom,len(startCodonPosList),len(stopCodonPosList))
    print (startCodonPosList[-10:])
    print (stopCodonPosList[-10:])
    
    findSorfsCoords(aChrom,aFrame,startCodonPosList,stopCodonPosList)










    