import sys
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--indexes', dest='inFile', required=True,
                    help='Path to indexes file, one per line.')
parser.add_argument('--config', dest='config', required=True)

o=parser.parse_args()

outFile = open(f"{o.config}.summary.csv",'w')

indexes = []
indexFile = open(o.inFile, 'r')
for line in indexFile:
    indexes.append(line.strip())
indexFile.close()
outFile.write("RunID,Index,Raw reads,SSCS reads,Mapped SSCS,DCS reads,Mapped DCS,% mapped SSCS,% mapped DCS,Raw/SSCS,SSCS/DCS,Peak Family Size,Max Family Size,Mean Insert Size,SSCS On Target,DCS On Target,DCS Mean Depth,DCS Max Depth,DCS Uncovered Target,Nucleotides Sequenced,A's Sequenced,T's Sequenced,C's Sequenced,G's Sequenced,Mutations,Mutation Frequency,A>T,A>C,A>G,T>A,T>C,T>G,C>A,C>T,C>G,G>A,G>T,G>C,ins,dels\n")

for index in indexes:
    print(f"Index {index}")
    print("Reading config")
    # Get the run ID from the config file
    runID=""
    c=""
    C=""
    d=""
    configFile = open(f"{index}/{index}_config.sh", 'r')
    for line in configFile:
        if "RUN_ID=" in line:
            runID = line.strip().split('=')[1].strip('"')
        elif "minClonal=" in line:
            c=line.strip().split('=')[1].split()[0]
        elif "maxClonal=" in line:
            C=line.strip().split('=')[1].split()[0]
        elif "minDepth=" in line:
            d=line.strip().split('=')[1].split()[0]
    configFile.close()
    
    print("Getting read counts")
    # get read counts
    # Read tagstats files:
    rawFlagstats = open(f"{index}/Stats/data/{runID}.temp.sort.flagstats.txt", 'r').readlines()
    sscsFlagstats = open(f"{index}/Stats/data/{runID}_mem.sscs.sort.flagstats.txt", 'r').readlines()
    dcsFlagstats = open(f"{index}/Stats/data/{runID}_mem.dcs.sort.flagstats.txt", 'r').readlines()
    rawReads = float(rawFlagstats[0].split()[0])
    #rawReads = float(pysam.flagstat(f"{index}/{runID}.temp.sort.bam").split('\n')[0].split()[0])
    #sscsFlagstat=pysam.flagstat(f"{index}/{runID}_mem.sscs.sort.bam").split('\n')
    sscsReads=float(sscsFlagstats[0].split()[0])
    mappedSscs=float(sscsFlagstats[4].split()[0])
    # ~ dcsFlagstat=pysam.flagstat(f"{index}/{runID}_mem.dcs.sort.bam").split('\n')
    dcsReads=float(dcsFlagstats[0].split()[0])
    mappedDcs=float(dcsFlagstats[4].split()[0])
    
    print("Processing Tagstats")
    # get tagstats numbers
    tagstatsFile = open(f"{index}/Stats/data/{runID}.tagstats.txt", 'r')
    lastProportion=1
    peakProportion = 0
    peakSize = 1
    maxSize=0
    for line in tagstatsFile:
        if float(line.split()[2]) <= lastProportion:
            lastProportion = float(line.split()[2])
        elif float(line.split()[2]) >= peakProportion:
            lastProportion = 0
            peakSize = line.split()[0]
            peakProportion = float(line.split()[2])
        maxSize = line.split()[0]
    tagstatsFile.close()
    sscsOnTarget="NA"
    
    # read depth file:
    print("Processing Depth")
    depthFile = open(f"{index}/Stats/data/{runID}.dcs.region.mutpos.vcf_depth.txt", 'r')
    totDepth = 0
    numLocs = 0
    dcsMaxDepth = 0
    for line in depthFile:
        if "#" not in line:
            totDepth += int(line.split('\t')[3])
            numLocs += 1
            dcsMaxDepth = max(dcsMaxDepth, int(line.split('\t')[3]))
    dcsOnTarget="NA"
    if numLocs != 0:
        dcsMeanDepth=totDepth / numLocs
    else:
        dcsMeanDepth=0
    dcsUncovered="NA"
    depthFile.close()
    # insert size file
    print("Processing Insert Size")
    insertSizeFile = open(f"{index}/Stats/data/{runID}.dcs.iSize_Metrics.txt", 'r')
    totInsertSize = 0
    numInsertReads = 0
    line = next(insertSizeFile)
    while "## HISTOGRAM" not in line:
        line = next(insertSizeFile)
    contIter = True
    line = next(insertSizeFile)
    while contIter:
        try:
            line = next(insertSizeFile)
            if line.strip() != "":
                linebins = [int(x) for x in line.strip().split('\t')]
                totInsertSize += linebins[0] * linebins[1]
                numInsertReads += linebins[1]
        except StopIteration:
            contIter = False
    if numInsertReads == 0:
        meanInsertSize = "N/A"
    else:
        meanInsertSize = totInsertSize / numInsertReads
    print("Processing countmuts")
    # get countmuts data
    sys.stderr.write(f"{index}/{runID}.dcs.filt.no_overlap.region.c{c}-{C}.d{d}.unique.countmuts.txt\n")

    cmFile = open(f"{index}/Final/dcs/{runID}.dcs.countmuts.csv", 'r')
    AsSeq=""
    AtoT=""
    AtoC=""
    AtoG=""
    TsSeq=""
    TtoA=""
    TtoC=""
    TtoG=""
    CsSeq=""
    CtoA=""
    CtoT=""
    CtoG=""
    GsSeq=""
    GtoA=""
    GtoT=""
    GtoC=""
    totalNt=""
    totalMuts=""
    ins=""
    dels=""
    
    for line in cmFile:
        if "##" not in line and "OVERALL" in line:
            linebins = line.strip().split(',')
            if "A>T" in line:
                AtoT=linebins[4]
                AsSeq=linebins[5]
            elif "A>C" in line:
                AtoC=linebins[4]
            elif "A>G" in line:
                AtoG=linebins[4]
            elif "T>A" in line:
                TtoA=linebins[4]
                TsSeq=linebins[5]
            elif "T>C" in line:
                TtoC=linebins[4]
            elif "T>G" in line:
                TtoG=linebins[4]
            elif "C>A" in line:
                CtoA=linebins[4]
                CsSeq=linebins[5]
            elif "C>T" in line:
                CtoT=linebins[4]
            elif "C>G" in line:
                CtoG=linebins[4]
            elif "G>A" in line:
                GtoA=linebins[4]
                GsSeq=linebins[5]
            elif "G>T" in line:
                GtoT=linebins[4]
            elif "G>C" in line:
                GtoC=linebins[4]
            elif "Total" in line and "SNV" in line:
                totalNt = float(linebins[5])
                totalMuts = float(linebins[4])
            elif "Total" in line and "INS" in line:
                ins=linebins[4]
            elif "Total" in line and "DEL" in line:
                dels=linebins[4]
            
    cmFile.close()
    if sscsReads > 0:
        percentMappedSSCS = mappedSscs/sscsReads
        rawPerSSCS = rawReads/sscsReads
    else:
        percentMappedSSCS = 0
        rawPerSSCS = 0
    if dcsReads > 0:
        percentMappedDCS = mappedDcs/dcsReads
        sscsPerDCS = sscsReads/dcsReads
    else:
        percentMappedDCS = 0
        sscsPerDCS = 0
    if totalNt > 0:
        mutFreq = totalMuts/totalNt
    else:
        mutFreq = 0
    outFile.write(
        f"{runID},"
        f"{index},{rawReads},{sscsReads},{mappedSscs},{dcsReads},{mappedDcs},"
        f"{percentMappedSSCS},{percentMappedDCS},{rawPerSSCS},{sscsPerDCS},"
        f"{peakSize},{maxSize},{meanInsertSize},{sscsOnTarget},{dcsOnTarget},{dcsMeanDepth},"
        f"{dcsMaxDepth},{dcsUncovered},{totalNt},{AsSeq},{TsSeq},{CsSeq},{GsSeq},{totalMuts},{mutFreq},"
        f"{AtoT},{AtoC},{AtoG},{TtoA},{TtoC},{TtoG},{CtoA},{CtoT},{CtoG},{GtoA},"
        f"{GtoT},{GtoC},{ins},{dels}\n"
        )
        
outFile.close()
        
