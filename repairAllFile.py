import glob, os
for bedFileName in glob.glob("output/*.bed"):
    print(bedFileName)
    bedFileNameOld=bedFileName+".old"
    os.rename(bedFileName,bedFileNameOld)
    bedFileOld=open(bedFileNameOld,"r")
    bedFile=open(bedFileName,"w")
    for aLine in bedFileOld:
        if not aLine[:3]==">>>":
            # print (aLine)
            bedFile.write(aLine)
bedFileOld.close()
bedFile.close()




