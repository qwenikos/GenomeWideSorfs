import sys
from  genomeScan_lib import *
from encoding_lib import *

########################################################################
########################################################################
########################################################################
########################################################################


###Global Vars
codonList = [
    "TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG",
    "ATT", "ATC", "ATA", "ATG", "GTT", "GTC", "GTA", "GTG",
    "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG",
    "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG",
    "TAT", "TAC", "TAA", "TAG", "CAT", "CAC", "CAA", "CAG",
    "AAT", "AAC", "AAA", "AAG", "GAT", "GAC", "GAA", "GAG",
    "TGT", "TGC", "TGA", "TGG", "CGT", "CGC", "CGA", "CGG",
    "AGT", "AGC", "AGA", "AGG", "GGT", "GGC", "GGA", "GGG"
]
# codonList = ["ATG"]

##Stop condon are three TAG TAA TGA and some alternatives (to check if exist for human)
allStartCodonList=["AAG","ACG","AGG","ATA","ATC","ATG","ATT","CTG","GTG","TTG"]
canonicalStartCodonList=["ATG"] ##canonical
nonCanonicalStartCodonList=["AAG","ACG","AGG","ATA","ATC","ATT","CTG","GTG","TTG"] ##non-canonical
stopCodonList=["TAG","TAA","TGA"]
# stopCodonList=["TGA"]

chromList=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
chromList=["22"]

frameList=["1","2","3","-1","-2","-3"]
frameList=["1"]

dictFileName="data/dnaDict.pkl"

num_processes=15
####################### FUNCTIONS ################################################

dnaDict=readDict(dictFileName)

def myFunction(aChrom,dnaSeq,codonList,frame):
    print(f"called for chrom={aChrom} and codon={codonList} and frame={frame}")

    dnaSeq=getSeqFrame(dnaSeq,frame)
    dnaSeqLen=len(dnaSeq)
    codonPosList, codonPosListLen  = findCodonPos(dnaSeq,codonList)
    # codonListStr=",".join(codonList)
    codon=codonList[0]
    aa=codonDNA_to_amino_acid[codon]
    outputLineList=[aChrom,codon,aa,dnaSeqLen,frame,codonPosListLen]
    outputLine="\t".join(map(str,outputLineList))
    return(outputLine)

########################### MAIN ######################################
headerList=["chromName","chromLength","frame","numbersOfATG","numberOfNonCanonicalCodons","numberOfStopCodons"]

test=False

seq=readChromFromFastaToSeq("data/Ensembl/Homo_sapiens.GRCh38.dna.chromosome.22.fa")

### create args for runs
argsList=[]
for aChrom in chromList:
    chromDnaSeq=dnaDict[aChrom].upper()
    # chromDnaSeq=seq
    print ("working on chromosome",aChrom)
    for aFrame in frameList:
        for aCodon in codonList:
            selectedCodonList=[aCodon] 
            argsTuple=(aChrom,chromDnaSeq,selectedCodonList,aFrame)
            argsList+=[argsTuple]

timeStart()

results=run_function_multiple_times(myFunction, argsList, num_processes)
timeEnd()

####### write results to file
outputFileName="allCodonStats.csv"
outputFile=open(outputFileName,"w")
headerLine="\t".join(map(str,headerList))
print (headerLine)
outputFile.write(headerLine+"\n")

for i, result in enumerate(results):
    print(result)
    outputFile.write(result+"\n")

outputFile.close()
