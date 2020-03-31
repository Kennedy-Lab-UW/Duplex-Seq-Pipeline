# BamToCountMuts.py
import pysam
from collections import Counter, defaultdict
from argparse import ArgumentParser
import sys
import re
import csv
import string
import logging
import gzip
import datetime
from math import sqrt

def Wilson(positive,  total) :

        if total == 0:

            return 0

        z = 1.96 #1.96 = 95%
        phat = float(positive) / total
        positiveCI = (phat 
                      + z * z / (2 * total) 
                      + z * sqrt(
                          (phat * (1 - phat) + z * z / (4 * total)) / total
                          )
                      ) / ( 1 + z * z / total)
        negativeCI =  (phat 
                       + z * z / (2 * total) 
                       - z * sqrt(
                           (phat * (1 - phat) + z * z / (4 * total)) / total
                           )
                       ) / (1 + z * z / total)

        return  (phat, positiveCI , negativeCI )

class Bed_File:
    def __init__(self, inFile):
        self.file = open(inFile, 'r')
        
    def __iter__(self):
        return(self)
    
    def __next__(self):
        myLine = self.file.readline()
        if myLine != "":
            return(Bed_Line(*myLine.strip().split()))
        else:
            raise StopIteration
            
class Bed_Line:
    def __init__(self, 
                 chrom, 
                 startPos, 
                 endPos, 
                 name="",
                 score="", 
                 strand="+",
                 thickStart="", 
                 thickEnd="",
                 itemRGB="", 
                 blockCount="",
                 blockSizes="", 
                 blockStarts=""):
        self.chrom = chrom
        self.startPos = int(startPos)
        self.endPos = int(endPos)
        if name == "":
            self.name = f"{self.chrom}:{self.startPos}-{self.endPos}"
        else:
            self.name = name
        if score == '':
            self.score = 0
        else:
            self.score = int(score)
        if strand in (".","+","-"):
            self.strand = strand
        else:
            self.strand = "."
        if thickStart == "":
            self.thickStart = self.startPos
        else:
            self.thickStart = int(thickStart)
        if thickEnd == "":
            self.thickEnd = self.endPos
        else:
            self.thickEnd = int(thickEnd)
        if itemRGB == "":
            self.itemRBG = None
        else:
            self.itemRGB = itemRGB
        if blockCount == "":
            self.blockCount = 1
        else:
            self.blockCount = int(blockCount)
        if blockSizes == "":
            self.blockSizes = [abs(self.endPos - self.startPos)]
        else:
            self.blockSizes = [int(x) for x in blockSizes.split(',')]
        if blockStarts == "":
            self.blockStarts = [0]
        else:
            self.blockStarts = [int(x) for x in blockStarts.split(',')]
    
    def samtoolsStr(self):
        return(f"{self.chrom}:{self.startPos}:{self.endPos}")
    
    def get_subregions(self):
        if self.strand == '+':
            return([Bed_Line(
                self.chrom, 
                self.startPos + self.blockStarts[x], 
                self.startPos + self.blockStarts[x] + self.blockSizes[x], 
                f"{self.name}_block{x + 1}", 
                self.score, 
                self.strand
                ) for x in range(len(self.blockStarts))])
        elif self.strand == '-':
            return([Bed_Line(
                self.chrom, 
                self.startPos + self.blockStarts[-x - 1], 
                self.startPos + self.blockStarts[-x - 1] + self.blockSizes[-x - 1], 
                f"{self.name}_block{x + 1}", 
                self.score, 
                self.strand
                ) for x in range(len(self.blockStarts))])
        
    def contains(self, inChr, inPos):
        if (
                inChr == self.chrom
                and inPos >= self.startPos
                and inPos < self.endPos
                ):
            return(True)
        else:
            return(False)
    
def getParams():
    parser = ArgumentParser(
        description=(
            f"Make a mutpos vcf file from a post-procesing duplex "
            f"sequencing BAM file.  Note that this program will fail if"
            f" the maximum depth in your bam file is > 8000.  "
            )
        )
    parser.add_argument(
        '-i', '--inBam', 
        action ='store', 
        dest = 'inBam', 
        help = 'An imput bam file. If None, defaults to stdin. [%(default)s]', 
        default = None
        )
    parser.add_argument(
        '-b', '--inBed', 
        action ='store', 
        dest = 'inBed', 
        help = 'An input bed file. If None, processes all positions. [%(default)s]', 
        default = None
        )
    parser.add_argument(
        '-f', '--inFasta', 
        action ='store', 
        dest = 'in_fasta', 
        help = 'An input bed file. If None, processes all positions. [%(default)s]', 
        default = None
        )
    parser.add_argument(
        '-o', '--outfile', 
        action = 'store', 
        dest = 'out_file', 
        help = 'A filename for the output file.  If None, outputs to stdout.  [%(default)s]', 
        default = None
        )
    parser.add_argument(
        "--round", 
        action="store", 
        type=int, 
        dest="round", 
        help="How many digits to round frequencies to.", 
        default = 4
        )
    parser.add_argument(
        '--logLevel', 
        action="store", 
        dest="logLvl", 
        default="Info",
        help=(f"Identification for how much information gets output. "
              f"Acceptable levels are: 'DEBUG', 'INFO', 'WARNING', "
              f"'ERROR', and 'CRITICAL'.  "
              )
        )
    parser.add_argument(
        '--samp_name', 
        action='store', 
        dest='sampName',
        default=None,
        help='A name for the sample in the output VCF file.  '
        )
    parser.add_argument(
        "-d", "--depth", 
        action = "store", 
        type = int, 
        dest = "mindepth", 
        default = 20,
        help = "Minimum depth for counting mutations at a site [20]"
        )
    parser.add_argument(
        "-c", "--min_clonality", 
        action = "store", 
        type = float, 
        dest = "min_clonality", 
        default = 0,
        help = (f"Minimum (exclusive) cutoff of mutant reads for scoring a "
                f"clonal mutation [0]"
                )
        )
    parser.add_argument(
        "-C", "--max_clonality", 
        action = "store", 
        type = float, 
        dest = "max_clonality", 
        default = 0.1,
        help = (f"Maximum (inclusive) cutoff of mutant reads for scoring a "
                f"clonal mutation [%(default)s]"
                )
        )
    parser.add_argument(
        "-n", "--n_cutoff", 
        action = "store", 
        type = float, 
        dest = "n_cutoff", 
        default = 0.05,
        help = "Maximum fraction of N's allowed to score a position [0.05]"
        )
    parser.add_argument(
        '-u', '--unique', 
        action = 'store_true', 
        dest = 'unique', 
        help = 'Run countMutsUnique instead of countMuts'
        )
    return(parser.parse_args())

class countMutsEngine:
    def __init__(self, 
                 inBam, 
                 inFasta, 
                 inBed=None, 
                 unique=False, 
                 Nprop=1, 
                 minDepth=1, 
                 minC=0, 
                 maxC=1, 
                 sampName=None
                 ):
        self.params = {
           "inBam":inBam, 
           "inFasta":inFasta, 
           "inBed":inBed, 
           "unique": unique, 
           "Nprop": Nprop, 
           "minDepth": minDepth, 
           "minC": minC, 
           "maxC": maxC
           }
        self.inBam = pysam.AlignmentFile(inBam, "rb")
        if sampName is None:
            self.sample = "Sample"
        else:
            self.sample = sampName
        self.minC = minC
        self.maxC = maxC
        self.inFasta = pysam.FastaFile(inFasta)
        self.rmNumTable = {ord(c): None for c in '1234567890'}
        if inBed is None:
            self.myBed = None
        else:
            self.myBed = Bed_File(inBed)
        self.mutsCounts = {
            "Aseq": 0, 
            "A>T": 0,
            "A>C": 0,
            "A>G": 0,
            "Tseq": 0,
            "T>A": 0,
            "T>C": 0,
            "T>G": 0,
            "Cseq": 0,
            "C>A": 0,
            "C>T": 0,
            "C>G": 0,
            "Gseq": 0,
            "G>A": 0,
            "G>T": 0,
            "G>C": 0,
            "ins":defaultdict(int),
            "dels": defaultdict(int),
            "DP": 0
            }
        self.subregCounts = {}
    
    def processLines(self, 
                     roundLevel
                     ):
        linesProcessed = 0
        if self.myBed is None:
            for pileup_column in self.inBam.pileup(
                    fastafile=self.inFasta, 
                    truncate=True, 
                    stepper="nofilter", 
                    max_depth=1000000, 
                    min_base_quality = 0
                    ):
                
                lnCnts = self.CountLine(
                    pileup_column, 
                    roundLevel,
                    self.params["unique"],
                    self.params["Nprop"], 
                    self.params["minDepth"], 
                    self.params["minC"], 
                    self.params["maxC"]
                    )
                self.mutsCounts["DP"] += lnCnts["DP"]
                self.mutsCounts[f"{lnCnts['RefBase']}seq"] += lnCnts["DP"]
                # ~ print(f"{pileup_column.reference_pos + 1}\t{lnCnts['DP']}\t"
                      # ~ f"{lnCnts['A'] + lnCnts['G'] + lnCnts['C'] + lnCnts['T']}")
                for xIter in ("A","C","G","T"):
                    if xIter != lnCnts["RefBase"]:
                        self.mutsCounts[f"{lnCnts['RefBase']}>{xIter}"] += lnCnts[xIter]
                for xIter in lnCnts["ins"]:
                    self.mutsCounts["ins"][xIter] += lnCnts["ins"][xIter]
                for xIter in lnCnts["dels"]:
                    self.mutsCounts["dels"][xIter] += lnCnts["dels"][xIter]
                
            linesProcessed += 1
            if linesProcessed % 1000 == 0:
                logging.info(f"{linesProcessed} lines processed...")
        else:
            linesCounted = []
            for myRegion in self.myBed:
                regStr = myRegion.samtoolsStr()
                subregions = myRegion.get_subregions()
                self.subregCounts[regStr] = {
                    "name": myRegion.name, 
                    "str": regStr, 
                    "Aseq": 0, 
                    "A>T": 0,
                    "A>C": 0,
                    "A>G": 0,
                    "Tseq": 0,
                    "T>A": 0,
                    "T>C": 0,
                    "T>G": 0,
                    "Cseq": 0,
                    "C>A": 0,
                    "C>T": 0,
                    "C>G": 0,
                    "Gseq": 0,
                    "G>A": 0,
                    "G>T": 0,
                    "G>C": 0,
                    "ins":defaultdict(int),
                    "dels": defaultdict(int),
                    "DP": 0
                    }
                if not ( 
                        len(subregions) == 1
                        and subregions[0].startPos == myRegion.startPos
                        and subregions[0].endPos == myRegion.endPos
                        ):
                    for subregion in subregions:
                        self.subregCounts[subregion.samtoolsStr()] = {
                            "name": subregion.name, 
                            "str": subregion.samtoolsStr(), 
                            "Aseq": 0, 
                            "A>T": 0,
                            "A>C": 0,
                            "A>G": 0,
                            "Tseq": 0,
                            "T>A": 0,
                            "T>C": 0,
                            "T>G": 0,
                            "Cseq": 0,
                            "C>A": 0,
                            "C>T": 0,
                            "C>G": 0,
                            "Gseq": 0,
                            "G>A": 0,
                            "G>T": 0,
                            "G>C": 0,
                            "ins":defaultdict(int),
                            "dels": defaultdict(int),
                            "DP": 0
                            }
                
                for pileup_column in self.inBam.pileup(
                        reference=myRegion.chrom, 
                        start=myRegion.startPos , 
                        end=myRegion.endPos,
                        fastafile=self.inFasta, 
                        truncate=True, 
                        stepper="nofilter", 
                        max_depth=1000000, 
                        min_base_quality = 0
                        ):
                    myChrPos = f"{myRegion.chrom}:{pileup_column.reference_pos + 1}"
                    lnCnts = self.CountLine(
                        pileup_column, 
                        roundLevel,
                        self.params["unique"],
                        self.params["Nprop"], 
                        self.params["minDepth"], 
                        self.params["minC"], 
                        self.params["maxC"]
                        )
                    self.subregCounts[regStr]["DP"] += lnCnts["DP"]
                    self.subregCounts[regStr][f"{lnCnts['RefBase']}seq"] += lnCnts["DP"]
                    
                    for xIter in ("A","C","G","T"):
                        if xIter != lnCnts["RefBase"]:
                            self.subregCounts[regStr][f"{lnCnts['RefBase']}>{xIter}"] += lnCnts[xIter]
                    for xIter in lnCnts["ins"]:
                        self.subregCounts[regStr]["ins"][xIter] += lnCnts["ins"][xIter]
                    for xIter in lnCnts["dels"]:
                        self.subregCounts[regStr]["dels"][xIter] += lnCnts["dels"][xIter]
                    
                    if not ( 
                            len(subregions) == 1
                            and subregions[0].startPos == myRegion.startPos
                            and subregions[0].endPos == myRegion.endPos
                            ):
                        for subregion in subregions:
                            if subregion.contains(myRegion.chrom, pileup_column.reference_pos):
                                self.subregCounts[subregion.samtoolsStr()]["DP"] += lnCnts["DP"]
                                self.subregCounts[subregion.samtoolsStr()][f"{lnCnts['RefBase']}seq"] += lnCnts["DP"]
                                
                                for xIter in ("A","C","G","T"):
                                    if xIter != lnCnts["RefBase"]:
                                        self.subregCounts[subregion.samtoolsStr()][f"{lnCnts['RefBase']}>{xIter}"] += lnCnts[xIter]
                                for xIter in lnCnts["ins"]:
                                    self.subregCounts[subregion.samtoolsStr()]["ins"][xIter] += lnCnts["ins"][xIter]
                                for xIter in lnCnts["dels"]:
                                    self.subregCounts[subregion.samtoolsStr()]["dels"][xIter] += lnCnts["dels"][xIter]
                    if myChrPos not in linesCounted:
                        linesCounted.append(myChrPos)
                        
                        # ~ print(f"{pileup_column.reference_pos + 1}\t{lnCnts['DP']}\t"
                          # ~ f"{lnCnts['A'] + lnCnts['G'] + lnCnts['C'] + lnCnts['T']}")
                        self.mutsCounts["DP"] += lnCnts["DP"]
                        self.mutsCounts[f"{lnCnts['RefBase']}seq"] += lnCnts["DP"]
                        
                        for xIter in ("A","C","G","T"):
                            if xIter != lnCnts["RefBase"]:
                                self.mutsCounts[f"{lnCnts['RefBase']}>{xIter}"] += lnCnts[xIter]
                        for xIter in lnCnts["ins"]:
                            self.mutsCounts["ins"][xIter] += lnCnts["ins"][xIter]
                        for xIter in lnCnts["dels"]:
                            self.mutsCounts["dels"][xIter] += lnCnts["dels"][xIter]
                        
                        linesProcessed += 1
                        if linesProcessed % 1000 == 0:
                            logging.info(f"{linesProcessed} lines processed...")
                    
    
    def CountLine(self, 
                  pileup_column, 
                  roundLevel = 4, 
                  unique=False, 
                  Nprop=1, 
                  minDepth=1, 
                  minC=0, 
                  maxC=1
                  ):
        myReads = Counter(
            [x.upper() 
                for x in pileup_column.get_query_sequences(
                    add_indels=True
                    )
                ]
            )
        
        mReads = defaultdict(int)

        # myReads should be something like
        myChrom = pileup_column.reference_name
        myPos = pileup_column.reference_pos + 1
        # I'll need to pull the reference base from the fasta file
        
        myRefBase = self.inFasta.fetch(reference=myChrom, start=myPos-1, end=myPos).upper()
        indelStrip = {ord(c): None for c in '1234567890'}
        logging.debug(f"Refernce is: {myRefBase}")
        myVcfLines = []
        myRefCount = myReads[myRefBase]
        myNCount = sum([myReads[x] for x in myReads if x[0] == "N"])
        myTotCount = sum([myReads[x] for x in myReads])
        clonalities = {x: myReads[x]/(myTotCount - myNCount) if myTotCount - myNCount > 0 else 0 for x in myReads}
        for x in myReads:
            if clonalities[x] >= minC and clonalities[x] <= maxC:
                if unique:
                    mReads[x] = 1
                else:
                    mReads[x] = myReads[x]
        if (
                myNCount / myTotCount <= Nprop and 
                myTotCount >= minDepth
                ):
            
            mutsDict = {"RefBase": myRefBase, 
                        "A": mReads["A"], 
                        "C": mReads["C"], 
                        "G": mReads["G"], 
                        "T": mReads["T"], 
                        "ins": defaultdict(int), 
                        "dels": defaultdict(int), 
                        "DP": myTotCount-myNCount
                        }
            for readTypeKey in mReads:
                # ~ if '-' in readTypeKey:
                    # ~ print(readTypeKey)
                if readTypeKey != "*":
                    if '-' in readTypeKey:
                        myRefBPs = self.inFasta.fetch(reference=myChrom, start=myPos-1, end=myPos+len(readTypeKey[2:])).upper()
                        if 'N' not in myRefBPs:
                            mutsDict["dels"][''.join(c for c in readTypeKey[1:] if c.isdigit())] += mReads[readTypeKey]
                    elif '+' in readTypeKey:
                        if 'N' not in readTypeKey:
                            mutsDict["ins"][''.join(c for c in readTypeKey[1:] if c.isdigit())] += mReads[readTypeKey]
        else:
            mutsDict = {"RefBase": myRefBase, 
                        "A": 0, 
                        "C": 0, 
                        "G": 0, 
                        "T": 0, 
                        "ins": defaultdict(int), 
                        "dels": defaultdict(int), 
                        "DP": 0
                        }
        #~ logging.debug(f"MyLine is {str(myLine)}")
        return(mutsDict)
    
    def genSummary(self, Fout):
        logging.debug("Generating summary")
        logging.debug(self.mutsCounts)
        logging.debug(self.subregCounts)
        outFile = open(Fout, 'w')
        outFile.write(
            f"##CountMuts output:\n"
            f"##Input file: \t{self.params['inBam']}\n"
            f"##Input reference:\t{self.params['inFasta']}\n"
            f"##Input bed:\t{self.params['inBed']}\n"
            f"##Minimum Depth: \t{self.params['minDepth']}\n"
            f"##Clonality: \t{self.params['minC']}-{self.params['maxC']}\n"
            )
        if self.params["unique"]:
            outFile.write(
                "##Unique mutations only\n"
                )
        outFile.write("#SAMPLE,REGION,MUTATION_TYPE,MUTATION_CLASS,COUNT,DENOMINATOR,FREQUENCY\n")
        # ~ outFile.write(
            # ~ f"\n##OVERALL: \n"
            # ~ )
        # Overall Output:
        totPointMuts = 0

        for x in ("A","T","C","G"):
            # ~ outFile.write(
                # ~ f"\n{x}'s sequenced:\t{self.mutsCounts[''.join((x,'seq'))]}\n"
                # ~ f"Mutation type\t#\tFrequency\t95% positive CI\t95% negative CI\n"
                # ~ )
            for y in ("A","T","C","G"):
                if x != y:
                    wilsonCI = Wilson(
                        self.mutsCounts[f"{x}>{y}"],  
                        max(self.mutsCounts[f"{x}seq"], 1)
                        )
                    totPointMuts += self.mutsCounts[f'{x}>{y}']
                    outFile.write(
                        f"{self.sample},"
                        f"OVERALL,{x}>{y},SNV,"
                        f"{self.mutsCounts[f'{x}>{y}']},"
                        f"{self.mutsCounts[f'{x}seq']},"
                        f"{wilsonCI[0]:.2e}\n"
                        )
        
        wilsonCI = Wilson(totPointMuts, max(self.mutsCounts['DP'], 1))
        
        outFile.write(
            f"{self.sample},"
            f"OVERALL,Total,SNV,"
            f"{totPointMuts},"
            f"{self.mutsCounts['DP']},"
            f"{wilsonCI[0]:.2e}\n"
            )
        
        # insertions:
        if self.mutsCounts['ins'] != {}:
            #print(self.mutsCounts['ins'])
            insKeys = sorted(int(x) for x in self.mutsCounts['ins'])
            # ~ outFile.write(
                # ~ "\nMutation type\t#\tFrequency\t95% positive CI\t95% negative CI\n"
                # ~ )
            for n in insKeys:
                if self.mutsCounts['ins'][str(n)] != 0:
                    wilsonCI = Wilson(
                        self.mutsCounts['ins'][str(n)], 
                        max(self.mutsCounts['DP'], 1)
                        )
                    outFile.write(
                        f"{self.sample},"
                        f"OVERALL,+{n},INS,{self.mutsCounts['ins'][str(n)]},"
                        f"{self.mutsCounts['DP']},"
                        f"{wilsonCI[0]:.2e}\n"
                        )
            totIns = sum([self.mutsCounts['ins'][x] for x in self.mutsCounts['ins']])
            wilsonCI = Wilson(
                totIns, 
                max(self.mutsCounts['DP'], 1)
                )
            outFile.write(
                f"{self.sample},"
                f"OVERALL,Total,INS,{totIns},"
                f"{self.mutsCounts['DP']},"
                f"{wilsonCI[0]:.2e}\n"
                )
        if self.mutsCounts['dels'] != {}:
            #print(self.mutsCounts['dels'])
            delsKeys = sorted(int(x) for x in self.mutsCounts['dels'])
            # ~ outFile.write(
                # ~ "\nMutation type\t#\tFrequency\t95% positive CI\t95% negative CI\n"
                # ~ )
            for n in delsKeys:
                if self.mutsCounts['dels'][str(n)] != 0:
                    wilsonCI = Wilson(
                        self.mutsCounts['dels'][str(n)], 
                        max(self.mutsCounts['DP'], 1)
                        )
                    outFile.write(
                        f"{self.sample},"
                        f"OVERALL,-{n},DEL,{self.mutsCounts['dels'][str(n)]},"
                        f"{self.mutsCounts['DP']},"
                        f"{wilsonCI[0]:.2e}\n"
                        )
            totDels = sum([self.mutsCounts['dels'][x] for x in self.mutsCounts['dels']])
            wilsonCI = Wilson(
                totDels, 
                max(self.mutsCounts['DP'], 1)
                )
            outFile.write(
                f"{self.sample},"
                f"OVERALL,Total,DEL,{totDels},"
                f"{self.mutsCounts['DP']},"
                f"{wilsonCI[0]:.2e}\n"
                )
                
        # Detail counts:
        for subreg in self.subregCounts:
            # ~ outFile.write(
                # ~ f"\n##Detail on {self.subregCounts[subreg]['name']} at {self.subregCounts[subreg]['str']}\n"
                # ~ )
            totPointMuts = 0
            for x in ("A","T","C","G"):
                # ~ outFile.write(
                    # ~ f"\n{x}'s sequenced:\t{self.subregCounts[subreg][''.join((x,'seq'))]}\n"
                    # ~ f"Mutation type\t#\tFrequency\t95% positive CI\t95% negative CI\n"
                    # ~ )
                for y in ("A","T","C","G"):
                    if x != y:
                        wilsonCI = Wilson(
                            self.subregCounts[subreg][f"{x}>{y}"],  
                            max(self.subregCounts[subreg][f"{x}seq"], 1)
                            )
                        totPointMuts += self.subregCounts[subreg][f'{x}>{y}']
                        outFile.write(
                            f"{self.sample},"
                            f"{self.subregCounts[subreg]['name']},"
                            f"{x}>{y},SNV,"
                            f"{self.subregCounts[subreg][f'{x}>{y}']},"
                            f"{self.subregCounts[subreg][f'{x}seq']},"
                            f"{wilsonCI[0]:.2e}\n"
                            )
            
            wilsonCI = Wilson(totPointMuts, max(self.subregCounts[subreg]['DP'], 1))
            outFile.write(
                f"{self.sample},"
                f"{self.subregCounts[subreg]['name']},"
                f"Total,SNV,"
                f"{totPointMuts},"
                f"{self.subregCounts[subreg]['DP']},"
                f"{wilsonCI[0]:.2e}\n"
                )
            
            # insertions:
            if self.subregCounts[subreg]['ins'] != {}:
                #print(self.subregCounts[subreg]['ins'])
                insKeys = sorted(int(x) for x in self.subregCounts[subreg]['ins'])
                # ~ outFile.write(
                    # ~ "\nMutation type\t#\tFrequency\t95% positive CI\t95% negative CI\n"
                    # ~ )
                for n in insKeys:
                    if self.subregCounts[subreg]['ins'][str(n)] != 0:
                        wilsonCI = Wilson(
                            self.subregCounts[subreg]['ins'][str(n)], 
                            max(self.subregCounts[subreg]['DP'], 1)
                            )
                        outFile.write(
                            f"{self.sample},"
                            f"{self.subregCounts[subreg]['name']},"
                            f"+{n},INS,{self.subregCounts[subreg]['ins'][str(n)]},"
                            f"{self.subregCounts[subreg]['DP']},"
                            f"{wilsonCI[0]:.2e}\n"
                            )
                totIns = sum([self.subregCounts[subreg]['ins'][x] for x in self.subregCounts[subreg]['ins']])
                wilsonCI = Wilson(
                    totIns, 
                    max(self.subregCounts[subreg]['DP'], 1)
                    )
                outFile.write(
                    f"{self.sample},"
                    f"{self.subregCounts[subreg]['name']},"
                    f"Total,INS,{totIns},"
                    f"{self.subregCounts[subreg]['DP']},"
                    f"{wilsonCI[0]:.2e}\n"
                    )
            if self.subregCounts[subreg]['dels'] != {}:
                #print(self.subregCounts[subreg]['dels'])
                delsKeys = sorted(int(x) for x in self.subregCounts[subreg]['dels'])
                # ~ outFile.write(
                    # ~ "\nMutation type\t#\tFrequency\t95% positive CI\t95% negative CI\n"
                    # ~ )
                for n in delsKeys:
                    if self.subregCounts[subreg]['dels'][str(n)] != 0:
                        wilsonCI = Wilson(
                            self.subregCounts[subreg]['dels'][str(n)], 
                            max(self.subregCounts[subreg]['DP'], 1)
                            )
                        outFile.write(
                            f"{self.sample},"
                            f"{self.subregCounts[subreg]['name']},"
                            f"+{n},DEL,{self.subregCounts[subreg]['dels'][str(n)]},"
                            f"{self.subregCounts[subreg]['DP']},"
                            f"{wilsonCI[0]:.2e}\n"
                            )
                totDels = sum([self.subregCounts[subreg]['dels'][x] for x in self.subregCounts[subreg]['dels']])
                wilsonCI = Wilson(
                    totDels, 
                    max(self.subregCounts[subreg]['DP'], 1)
                    )
                outFile.write(
                    f"{self.sample},"
                    f"{self.subregCounts[subreg]['name']},"
                    f"Total,DEL,{totDels},"
                    f"{self.subregCounts[subreg]['DP']},"
                    f"{wilsonCI[0]:.2e}\n"
                    )
                    
        outFile.close()

def main():
    o = getParams()
    numeric_level = getattr(logging, o.logLvl.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: {o.logLvl}')
    logging.basicConfig(
        format='%(levelname)s: %(message)s', 
        level=numeric_level, 
        )
    # Start engine
    logging.info("Starting Engine")
    myEngine = countMutsEngine(
        o.inBam, 
        o.in_fasta, 
        o.inBed, 
        o.unique,
        o.n_cutoff, 
        o.mindepth, 
        o.min_clonality, 
        o.max_clonality, 
        o.sampName
        )
    logging.info("Processing Lines")
    '''
        def processLines(self, 
                     roundLevel, 
                     unique=False, 
                     Nprop=1, 
                     minDepth=1, 
                     minC=0, 
                     maxC=1
    '''
    myEngine.processLines(
        o.round
        )
    myEngine.genSummary(o.out_file)
            
            
if __name__ == "__main__":
    main()
