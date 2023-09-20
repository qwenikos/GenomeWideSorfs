############################################
import time
start_time = time.time()
############################################

import pickle

# genomeFastaFileName="/home/nikos/Bioinfo/sorfs/data/GCF_000001405.40_GRCh38.p14_genomic.fna" ##NCBI
genomeFastaFileName="/home/nikos/Bioinfo/sorfs/data/Ensembl/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa"  ##Ensembl
genomeFastaFile=open(genomeFastaFileName,"r")
genomeFastaOneLineFileName="/home/nikos/Bioinfo/sorfs/data/GRCh38_oneLine.fasta"
genomeFastaOneLineFile=open(genomeFastaOneLineFileName,"w")

headerLine=""
seqStr=""
cnt=0
dnaDict={}

for aLine in genomeFastaFile:
    if aLine[0]==">":
        if seqStr=="": ##means first loop
            fistLoop=True
        else: 
            seqLen=len(seqStr)
            print (chrom,"\t",seqLen)
            if  len(chrom)<3:
                genomeFastaOneLineFile.write(headerLine+"|len="+str(seqLen)+"\n")
                genomeFastaOneLineFile.write(seqStr+"\n")
                dnaDict[chrom]=seqStr

        headerLine=aLine
        chrom=headerLine.rstrip().split(" ")[0][1:]
        seqStr=""
        
    else:
        seqStr=seqStr+aLine.rstrip()

print ("num of seq=",cnt)

genomeFastaFile.close()
genomeFastaOneLineFile.close()

serializedDnaDictFileName="data/dnaDict.pkl"
serializedDnaDict = pickle.dumps(dnaDict)
serializedDnaDictFile = open(serializedDnaDictFileName,"wb")
pickle.dump(dnaDict,serializedDnaDictFile)
serializedDnaDictFile.close()


################################################
print("--- %s seconds ---" % (time.time() - start_time))
################################################

