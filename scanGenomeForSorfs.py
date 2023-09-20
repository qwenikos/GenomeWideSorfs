import time
def timeStart():
    import time
    global start_time 
    start_time  = time.time()
    print ("++++++time start+++++")

#######################################################################

def timeEnd():
    print("--- %s seconds ---" % (time.time() - start_time))

#######################################################################

def readDict(fileName):
    import pickle
    serializedDnaDictFileName=fileName
    serializedDnaDictFile=open(serializedDnaDictFileName,"rb")
    dnaDict=pickle.load(serializedDnaDictFile)
    serializedDnaDictFile.close()
    return dnaDict

#######################################################################

def reverseComplement(dna):
    # print("---------------------")
    # print ("dna=",dna)
    # print("---------------------")
    dnaUpper = dna.upper()
    # complement strand
    # print (dnaUpper)
    dnaRev=""
    for nt in dna:
        if nt=="A":
            dnaRev+="T"
        elif nt=="T":
            dnaRev+="A"
        elif nt=="C":
            dnaRev+="G"
        elif nt=="G":
            dnaRev+="C"
        else:
            dnaRev+=nt
    # reverse strand
    dnaRevComp=dnaRev[::-1]
    # print("---------------------")
    # print(dnaRevComp)
    # print("---------------------")
    return dnaRevComp

def runFunction(aChrom,dna,frame):
    dnaSeq=""

    startCodonPosList=[]
    stopCodonPosList=[]
    # print ("dna=",dna)

    if frame=="1":
            dnaSeq=dna[0:]
    if frame=="2":
            dnaSeq=dna[1:]
    if frame=="3":
            dnaSeq=dna[2:]

    
    # print ("revCompDna=",revCompDna)
    # print ("revCompDna",revCompDna)
    if (frame=="-1"):
        revCompDna=reverseComplement(dna)
        dnaSeq=revCompDna[0:]
    if (frame=="-2"):
        revCompDna=reverseComplement(dna)
        dnaSeq=revCompDna[1:]
    if (frame=="-3"):
        revCompDna=reverseComplement(dna)
        dnaSeq=revCompDna[2:]
    # print (dnaSeq,"dnaSeq")
    print ("frame=",frame)
    
    dnaSeqlen=len(dnaSeq)
    # print (dnaSeq[-1000:])
    print ("dnaSeqlen",dnaSeqlen)
    for idx in range(0,dnaSeqlen,3):
        # print ("idx",idx)
        trinucleotide=dnaSeq[idx:idx+3] 
        # print (trinucleotide)
        if trinucleotide in startCodonList:
            startPos=idx #convert to coords
            startCodonPosList+=[startPos]
            # print (trinucleotide)
        if trinucleotide in sporCodonList:
            stopPos=idx+3 #convert to coords
            stopCodonPosList+=[stopPos]
    return startCodonPosList,stopCodonPosList

########################################################################

def findSorfsCoords(chrom,aFrame,startCodonPosList,stopCodonPosList):
    # startCodonPosList=tuple(startCodonPosList[:20000])
    # stopCodonPosList=tuple(stopCodonPosList[:20000])
    chromName="chr"+chrom
    if aFrame=="1":
        frameName="_plus_1"
    elif (aFrame=="-1"):
        frameName="_minus_1"
    elif aFrame=="2":
        frameName="_plus_2"
    elif (aFrame=="-2"):
        frameName="_minus_2"
    elif aFrame=="3":
        frameName="_plus_3"
    elif (aFrame=="-3"):
        frameName="_minus_3"
    else:
        print ("erron in frame")
        exit()
    outputFile=open("output/"+chromName+frameName+".bed","w")
    sorfCnt=0
    idxa=0
    idxb=0
    outputColsStr=""
    startCodonPosListLen=len(startCodonPosList)
    
    for aItem in startCodonPosList:

        idxa+=1
        idxPercent=round((idxa/startCodonPosListLen),4)*100
        if (idxPercent%20)==0:
            print (idxPercent,"%", "sorfs count=",sorfCnt)
            # outputFile.write(">>>>>>>>>>>>>>>>>"+str(idxPercent)+">>>>>>>>>>>>>>>\n")
            

        startPosForStops=0
        for bItem in stopCodonPosList[startPosForStops:]:
            # print (">>>>>>>>>>>>>>>>>")
            # print (bItem)
        # for bItem in stopCodonPosList:
            idxb+=1
            diff=(bItem-aItem)
            if (diff<0):
                continue
            elif (diff<=300) and (diff>=24):  ## At least 8AA
                startPosForStops=idxb ## update for improve perf
                sorfCnt+=1
                outputColsList=[chromName,str(aItem),str(bItem),"name","0","+"]
                outputColsStr="\t".join(outputColsList)+"\n"
                outputFile.write(outputColsStr) 
                # print ("sorf",aItem,bItem,diff, )
                break
            else:
                startPosForStops=idxb ## update for improve perf
                # print (">")
                break
       
    outputFile.close()
    # print (sorfCnt)

########################################################################
########################################################################
########################################################################
########################################################################


###Global Vars

##Stop condon are three TAG TAA TGA and some alternatives (to check if exist for human)
startCodonList=["AAG","ACG","AGG","ATA","ATC","ATG","ATT","CTG","GTG","TTG"]
startCodonList=["ATG"]

sporCodonList=["TAG","TAA","TGA"]
sporCodonList=["TGA"]

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










    