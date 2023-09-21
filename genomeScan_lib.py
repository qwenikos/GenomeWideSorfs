import time
import multiprocessing
#######################################################################

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
        if trinucleotide in stopCodonList:
            stopPos=idx+3 #convert to coords
            stopCodonPosList+=[stopPos]
    return startCodonPosList,stopCodonPosList

def getSeqFrame(dnaSeq,frame): ##frame get values (1,2,3,-1,-2,-3)
    frame=str(frame)
    if frame=="1":
        frameSeq=dnaSeq[0:]
    if frame=="2":
        frameSeq=dnaSeq[1:]
    if frame=="3":
        frameSeq=dnaSeq[2:]
    if (frame=="-1"):
        revCompDna=reverseComplement(dnaSeq)
        frameSeq=revCompDna[0:]
    if (frame=="-2"):
        revCompDna=reverseComplement(dnaSeq)
        frameSeq=revCompDna[1:]
    if (frame=="-3"):
        revCompDna=reverseComplement(dnaSeq)
        frameSeq=revCompDna[2:]
    # print (frameSeq,"frameSeq")
    # print ("frame=",frame)
    return (frameSeq)

########################################################################
def readChromFromFastaToSeq(fastaFileName):
    sequence=""
    fastaFile=open(fastaFileName,"r")
    header=fastaFile.readline()
    for aLine in fastaFile:
        sequence+=aLine.rstrip().upper()
    return sequence
        
########################################################################
def findCodonPos(dnaSeq,codonList):
 
    codonPosList=[]
    codonPosListLen=0
    dnaSeqlen=len(dnaSeq)

    for idx in range(0,dnaSeqlen,3):
        # print ("idx",idx)
        trinucleotide=dnaSeq[idx:idx+3] 
        # print (trinucleotide)
        if trinucleotide in codonList:
            startPos=idx #convert to coords
            codonPosList+=[startPos]
            # print (startPos,dnaSeq[startPos:startPos+20])
            codonPosListLen+=1

    return codonPosList,codonPosListLen

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
##############################################################################
### arguments as list of tuples eg arguments_list = [(1, 'A'), (2, 'B'), (3, 'C'), (4, 'D'), (5, 'E')] 

def run_function_multiple_times(func, args_list, num_processes):

    # Create a pool of processes
    pool = multiprocessing.Pool(processes=num_processes)

    # Use the pool to map the function to the argument list
    results=pool.starmap(func, args_list)

    # Close the pool to prevent any more tasks from being submitted
    pool.close()

    # Wait for all processes to complete
    pool.join()

    return results
