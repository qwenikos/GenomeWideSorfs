import sys
from  genomeScan_lib import *

########################################################################
########################################################################
########################################################################
########################################################################


###Global Vars

##Stop condon are three TAG TAA TGA and some alternatives (to check if exist for human)
allStartCodonList=["AAG","ACG","AGG","ATA","ATC","ATG","ATT","CTG","GTG","TTG"]
canonicalStartCodonList=["ATG"] ##canonical
nonCanonicalStartCodonList=["AAG","ACG","AGG","ATA","ATC","ATT","CTG","GTG","TTG"] ##non-canonical
stopCodonList=["TAG","TAA","TGA"]
# stopCodonList=["TGA"]

chromList=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
# chromList=["22"]

frameList=["1","2","3","-1","-2","-3"]

startCodonPosList=[]
stopCodonPosList=[]

dictFileName="data/dnaDict.pkl"

####################### FUNCTIONS ################################################

dnaDict=readDict(dictFileName)

def myFunction(aChrom,dnaSeq,codon,frame):
    c_StartCodonPosListLen=[]
    nc_StartCodonPosListLen=[]
    stopCodonPosListLen=[]

    dnaSeq=getSeqFrame(dnaSeq,frame)
    dnaSeqLen=len(dnaSeq)
    c_StartCodonPosList, c_StartCodonPosListLen     = findCodonPos(dnaSeq,canonicalStartCodonList)
    nc_StartCodonPosList, nc_StartCodonPosListLen   = findCodonPos(dnaSeq,nonCanonicalStartCodonList)
    stopCodonPosList, stopCodonPosListLen           = findCodonPos(dnaSeq,stopCodonList)

    outputLineList=[aChrom,dnaSeqLen,frame,c_StartCodonPosListLen,nc_StartCodonPosListLen,stopCodonPosListLen]
    outputLine="\t".join(map(str,outputLineList))
    return(outputLine)

########################### MAIN ######################################
headerList=["chromName","chromLength","frame","numbersOfATG","numberOfNonCanonicalCodons","numberOfStopCodons"]

test=False


if test:
    timeStart()
    num_processes=1
    dnaseq="AAAGTGATGAAAAAGATGTGA"
    results=run_function_multiple_times(myFunction, [("chrom",dnaseq,["ATG"],["AAG"],["TGA"],"2")], num_processes)
    timeEnd()
    print (results)
    exit()

### create args for runs
argsList=[]
for aChrom in chromList:
    chromDnaSeq=dnaDict[aChrom]
    print ("working on chromosome",aChrom)
    for aFrame in frameList:
        argsTuple=(aChrom,chromDnaSeq,canonicalStartCodonList,nonCanonicalStartCodonList,stopCodonList,aFrame)
        argsList+=[argsTuple]

timeStart()
num_processes=1
results=run_function_multiple_times(myFunction, argsList, num_processes)
timeEnd()

####### write results to file
outputFileName="codonStats.csv"
outputFile=open(outputFileName,"w")
headerLine="\t".join(map(str,headerList))
print (headerLine)
outputFile.write(headerLine+"\n")

for i, result in enumerate(results):
    print(result)
    outputFile.write(result+"\n")

outputFile.close()
