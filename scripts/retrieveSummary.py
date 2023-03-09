import sys
from argparse import ArgumentParser
import pandas as pd
from snakemake.utils import validate

parser = ArgumentParser()
parser.add_argument('--config', dest='config', required=True,
                    help='Path to DS config file')

o = parser.parse_args()

samples = pd.read_csv(o.config).set_index("sample", drop=False)
validate(samples, f"{sys.path[0]}/../DS_baseSchema.yaml")

def get_sample(sample):
    return samples.loc[sample, "sample"]

def get_rglb(sample):
    return samples.loc[sample, "rglb"]

def get_rgpl(sample):
    return samples.loc[sample, "rgpl"]

def get_rgpu(sample):
    return samples.loc[sample, "rgpu"]

def get_rgsm(sample):
    return samples.loc[sample, "rgsm"]

def get_reference(sample):
    return samples.loc[sample, "reference"]

def get_target_bed(sample):
    return samples.loc[sample, "target_bed"]

def get_blast_db(sample):
    return samples.loc[sample, "blast_db"]

def get_blast_db_path(sample):
    return f'{samples.loc[sample, "blast_db"]}.nal'

def get_target_taxon(sample):
    return samples.loc[sample, "targetTaxonId"]

def get_baseDir(sample):
    return samples.loc[sample, "baseDir"]

def get_in1(sample):
    return f'{samples.loc[sample, "baseDir"]}/{samples.loc[sample, "in1"]}'

def get_in2(sample):
    return f'{samples.loc[sample, "baseDir"]}/{samples.loc[sample, "in2"]}'

def get_mqFilt(sample):
    return samples.loc[sample, "mqFilt"]

def get_minMem(sample):
    return samples.loc[sample, "minMem"]

def get_maxMem(sample):
    return samples.loc[sample, "maxMem"]

def get_cutOff(sample):
    return samples.loc[sample, "cutOff"]

def get_nCutOff(sample):
    return samples.loc[sample, "nCutOff"]

def get_umiLen(sample):
    return samples.loc[sample, "umiLen"]

def get_spacerLen(sample):
    return samples.loc[sample, "spacerLen"]

def get_locLen(sample):
    return samples.loc[sample, "locLen"]

def get_readLen(sample):
    return samples.loc[sample, "readLen"]

def get_clipBegin(sample):
    return samples.loc[sample, "clipBegin"]

def get_clipEnd(sample):
    return samples.loc[sample, "clipEnd"]

def get_runSSCS(sample):
    return samples.loc[sample, "runSSCS"]

def get_runDCS(sample):
    return samples.loc[sample, "runDCS"]

def get_makeDCS(sample):
    return samples.loc[sample, "makeDCS"]

def get_cm_outputs(sample):
    return samples.loc[sample, "cm_outputs"]

def get_cm_sumTypes(sample):
    return samples.loc[sample, "cm_sumTypes"]

def get_minClonal(sample):
    return samples.loc[sample, "minClonal"]

def get_maxClonal(sample):
    return samples.loc[sample, "maxClonal"]

def get_minDepth(sample):
    return samples.loc[sample, "minDepth"]

def get_maxNs(sample):
    return samples.loc[sample, "maxNs"]

def get_cleanup(sample):
    return samples.loc[sample, "cleanup"]

def get_recovery(sample):
    return f'{sys.path[0]}/scripts/RecoveryScripts/{samples.loc[sample, "recovery"]}'


def main():
    outFile = open(f"{o.config}.summary.csv", 'w')

    outFile.write(
        "RunID,Index,Raw reads,"
        "SSCS reads,Mapped SSCS,"
        "DCS reads,Mapped DCS,"
        "% mapped SSCS,% mapped DCS,"
        "Raw/SSCS,SSCS/DCS,"
        "Peak Family Size,Max Family Size,Mean Insert Size,"
        "Raw On Target,SSCS On Target,DCS On Target,"
        "DCS Mean Depth,DCS Max Depth,"
        "DCS Uncovered Target,"
        "Nucleotides Sequenced,"
        "A's Sequenced,T's Sequenced,C's Sequenced,G's Sequenced,"
        "Mutations,Mutation Frequency,"
        "A>T,A>C,A>G,T>A,T>C,T>G,"
        "C>A,C>T,C>G,G>A,G>T,G>C,"
        "ins,dels\n"
    )

    for sampIter in samples.index:
        print(f"Sample {sampIter}")
        print("Reading config")
        # Get the run ID from the config file
        runID = sampIter
        c = get_minClonal(sampIter)
        C = get_minClonal(sampIter)
        d = get_minDepth(sampIter)
        baseDir = get_baseDir(sampIter)

        print("Getting read counts")
        # get read counts
        # Read tagstats files:
        rawFlagstats = open(f"{baseDir}/Stats/data/{runID}.temp.sort.flagstats.txt", 'r').readlines()
        sscsFlagstats = open(f"{baseDir}/Stats/data/{runID}_mem.sscs.sort.flagstats.txt", 'r').readlines()
        dcsFlagstats = open(f"{baseDir}/Stats/data/{runID}_mem.dcs.sort.flagstats.txt", 'r').readlines()
        rawReads = float(rawFlagstats[0].split()[0])
        # rawReads = float(pysam.flagstat(f"{index}/{runID}.temp.sort.bam").split('\n')[0].split()[0])
        # sscsFlagstat=pysam.flagstat(f"{index}/{runID}_mem.sscs.sort.bam").split('\n')
        sscsReads = float(sscsFlagstats[0].split()[0])
        mappedSscs = float(sscsFlagstats[6].split()[0])
        # ~ dcsFlagstat=pysam.flagstat(f"{index}/{runID}_mem.dcs.sort.bam").split('\n')
        dcsReads = float(dcsFlagstats[0].split()[0])
        mappedDcs = float(dcsFlagstats[6].split()[0])

        print("Processing Tagstats")
        # get tagstats numbers
        tagstatsFile = open(f"{baseDir}/Stats/data/{runID}.tagstats.txt", 'r')
        lastProportion = 1
        peakProportion = 0
        peakSize = 1
        maxSize = 0
        for line in tagstatsFile:
            if float(line.split()[2]) <= lastProportion:
                lastProportion = float(line.split()[2])
            elif float(line.split()[2]) >= peakProportion:
                lastProportion = 0
                peakSize = line.split()[0]
                peakProportion = float(line.split()[2])
            maxSize = line.split()[0]
        tagstatsFile.close()

        # read raw on target file
        rawTarget = open(f"{baseDir}/Stats/data/{runID}_onTargetCount.txt", 'r').readlines()
        if int(rawTarget[1].split()[0]) == 0:
            rawOnTarget=0
        else:
            rawOnTarget=f"{round(int(rawTarget[0].split()[0])/int(rawTarget[1].split()[0]),4) * 100}%"

        # read sscs on target file
        sscsTarget = open(f"{baseDir}/Stats/data/{runID}.sscs_onTargetCount.txt", 'r').readlines()
        if int(sscsTarget[1].split()[0]) == 0:
            sscsOnTarget=0
        else:
            sscsOnTarget=f"{round(int(sscsTarget[0].split()[0])/int(sscsTarget[1].split()[0]),4) * 100}%"
        # read dcs on target file
        dcsTarget = open(f"{baseDir}/Stats/data/{runID}.dcs_onTargetCount.txt", 'r').readlines()
        if int(dcsTarget[1].split()[0]) == 0:
            dcsOnTarget=0
        else:
            dcsOnTarget=f"{round(int(dcsTarget[0].split()[0])/int(dcsTarget[1].split()[0]),4) * 100}%"

        # read depth file:
        print("Processing Depth")
        depthFile = open(f"{baseDir}/Stats/data/{runID}.dcs.depth.txt", 'r')
        totDepth = 0
        numLocs = 0
        dcsMaxDepth = 0
        for line in depthFile:
            if "#" not in line:
                totDepth += int(line.split('\t')[3])
                numLocs += 1
                dcsMaxDepth = max(dcsMaxDepth, int(line.split('\t')[3]))
        if numLocs != 0:
            dcsMeanDepth = totDepth / numLocs
        else:
            dcsMeanDepth = 0
        dcsUncovered = "NA"
        depthFile.close()
        # insert size file
        print("Processing Insert Size")
        insertSizeFile = open(f"{baseDir}/Stats/data/{runID}.dcs.iSize_Metrics.txt", 'r')
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
        sys.stderr.write(f"{baseDir}/{runID}.dcs.filt.no_overlap.region.c{c}-{C}.d{d}.unique.countmuts.txt\n")

        cmFile = open(f"{baseDir}/Final/dcs/{runID}.dcs.countmuts.csv", 'r')
        AsSeq = ""
        AtoT = ""
        AtoC = ""
        AtoG = ""
        TsSeq = ""
        TtoA = ""
        TtoC = ""
        TtoG = ""
        CsSeq = ""
        CtoA = ""
        CtoT = ""
        CtoG = ""
        GsSeq = ""
        GtoA = ""
        GtoT = ""
        GtoC = ""
        totalNt = ""
        totalMuts = ""
        ins = ""
        dels = ""

        for line in cmFile:
            if "##" not in line and "OVERALL" in line:
                linebins = line.strip().split(',')
                if "A>T" in line:
                    AtoT = linebins[4]
                    AsSeq = linebins[5]
                elif "A>C" in line:
                    AtoC = linebins[4]
                elif "A>G" in line:
                    AtoG = linebins[4]
                elif "T>A" in line:
                    TtoA = linebins[4]
                    TsSeq = linebins[5]
                elif "T>C" in line:
                    TtoC = linebins[4]
                elif "T>G" in line:
                    TtoG = linebins[4]
                elif "C>A" in line:
                    CtoA = linebins[4]
                    CsSeq = linebins[5]
                elif "C>T" in line:
                    CtoT = linebins[4]
                elif "C>G" in line:
                    CtoG = linebins[4]
                elif "G>A" in line:
                    GtoA = linebins[4]
                    GsSeq = linebins[5]
                elif "G>T" in line:
                    GtoT = linebins[4]
                elif "G>C" in line:
                    GtoC = linebins[4]
                elif "Total" in line and "SNV" in line:
                    totalNt = float(linebins[5])
                    totalMuts = float(linebins[4])
                elif "Total" in line and "INS" in line:
                    ins = linebins[4]
                elif "Total" in line and "DEL" in line:
                    dels = linebins[4]

        cmFile.close()
        if sscsReads > 0:
            percentMappedSSCS = mappedSscs / sscsReads
            rawPerSSCS = rawReads / sscsReads
        else:
            percentMappedSSCS = 0
            rawPerSSCS = 0
        if dcsReads > 0:
            percentMappedDCS = mappedDcs / dcsReads
            sscsPerDCS = sscsReads / dcsReads
        else:
            percentMappedDCS = 0
            sscsPerDCS = 0
        if totalNt > 0:
            mutFreq = totalMuts / totalNt
        else:
            mutFreq = 0
        outFile.write(
            f"{runID},{baseDir},{rawReads},"
            f"{sscsReads},{mappedSscs},"
            f"{dcsReads},{mappedDcs},"
            f"{percentMappedSSCS},{percentMappedDCS},"
            f"{rawPerSSCS},{sscsPerDCS},"
            f"{peakSize},{maxSize},{meanInsertSize},"
            f"{rawOnTarget},{sscsOnTarget},{dcsOnTarget},"
            f"{dcsMeanDepth},{dcsMaxDepth},"
            f"{dcsUncovered},"
            f"{totalNt},"
            f"{AsSeq},{TsSeq},{CsSeq},{GsSeq},"
            f"{totalMuts},{mutFreq},"
            f"{AtoT},{AtoC},{AtoG},{TtoA},{TtoC},{TtoG},"
            f"{CtoA},{CtoT},{CtoG},{GtoA},{GtoT},{GtoC},"
            f"{ins},{dels}\n"
        )

    outFile.close()


if __name__ == "__main__":
    main()
