import sys
import math
from argparse import ArgumentParser
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
from VCF_Parser import *

refConvert = {"C": "C", "G": "C", "T": "T", "A": "T"}
compBase = {"C": "G", "G": "C", "T": "A", "A": "T", "N": "N"}


class MismatchCounter:
    def __init__(self, max_read_length, filters, filterSet, indelsSet, varsSet):
        self.mismatch_counts = [
            {"C>T": 0, "C>A": 0, "C>G": 0, "T>A": 0, "T>C": 0, "T>G": 0, "C>N": 0, "T>N": 0, "Count": 0}
            for i in range(max_read_length)
        ]
        self.indelCounts = 0
        self.filters = filters
        self.snps = filterSet
        self.indels = indelsSet
        self.vars = varsSet
        #self.debug_file = open("MPC_debug.txt",'w')

    def countRead(self, read):
        numMuts = 0
        if read.has_tag("MD"):
            refChrom = read.reference_name
            # Get read direction
            if read.is_reverse:
                readDir = -1
            else:
                readDir = 1
            # get actual read length
            rLen = read.infer_read_length()
            hardclipping = 0
            if read.cigartuples[0][0] == 5:
                hardclipping = read.cigartuples[0][1]
            if rLen > len(self.mismatch_counts):
                sys.stderr.write("ERROR: Actual read length longer than given read length! \n")
                raise Exception
            readMD = read.get_aligned_pairs(False, True)
            md2 = [x for x in readMD if x[0] is not None and x[2] is not None]
            for x in md2:
                if (
                        read.query_sequence[x[0]] != x[2].upper()
                ):
                    refBase = x[2].upper()
                    readBase = read.query_sequence[x[0]]
                    refPos = x[1] + 1

                    if refBase not in ("C", "T"):
                        refBase = compBase[refBase]
                        readBase = compBase[readBase]
                    snpTest = f"{refChrom}:{refPos}:{refBase}>{readBase}"
                    use_variant = True
                    
                    if snpTest in self.snps:
                        # check for and remove SNPs if requested
                        use_variant = False
                    if (
                            any([near_indel(snpTest, x) for x in self.indels])
                            and "INDEL" in self.filters):
                        # check for and remove variants near indels, if requested
                        use_variant = False
                    if snpTest not in self.vars and "INCLUDE" in self.filters:
                        # remove non-included variants, if requested
                        use_variant = False
                    if use_variant:
                        #self.debug_file.write(f"{snpTest}\n")
                        if readDir == -1:
                            self.mismatch_counts[rLen - x[0] - hardclipping - 1][f"{refBase}>{readBase}"] += 1
                        elif readDir == 1:
                            self.mismatch_counts[x[0] + hardclipping][f"{refBase}>{readBase}"] += 1
                        if readBase != "N":
                            numMuts += 1
                if readDir == -1:
                    self.mismatch_counts[rLen - x[0] - hardclipping - 1]["Count"] += 1
                elif readDir == 1:
                    self.mismatch_counts[x[0] + hardclipping]["Count"] += 1
            indelTest1 = [x[0] for x in read.cigartuples]
            if 1 in indelTest1 or 2 in indelTest1:
                indelCheck = getIndels(read)
                for indelIter in indelCheck:
                    if (
                            indelIter not in self.snps
                            and 'N' not in indelIter.split(':')[2].split('>')[0]
                            and 'N' not in indelIter.split(':')[2].split('>')[1]
                    ):
                        numMuts += 1
        return numMuts

    def getCounts(self):
        return self.mismatch_counts


def getIndels(read):
    contLoop = False
    if read.has_tag('MD'):
        md = read.get_aligned_pairs(False, True)
        readSeq = read.query_sequence
        contLoop = True
        x = 0
        # Get start softclipping
        while md[x][2] is None and contLoop:
            x += 1
            if x == len(md):
                contLoop = False
        # Process center
        contLoop = True
        indelDict = []
        readChr = read.reference_name
        while contLoop:
            indelSeq = [""]
            indelRef = ""
            indelPos = 0
            # test for insertion
            if md[x][2] is None:
                indelSeq.append(readSeq[md[x - 1][0]])
                indelRef = md[x - 1][2]
                indelPos = md[x - 1][1] + 1
                while contLoop and md[x][2] is None:
                    indelSeq.append(readSeq[md[x][0]])
                    x += 1
                    if x == len(md):
                        contLoop = False
                if contLoop:
                    indelDict.append(
                        f"{readChr}:"
                        f"{indelPos}:"
                        f"{indelRef}>{''.join(indelSeq)}"
                    )
            elif md[x][0] is None:
                indelSeq.append(md[x - 1][2])
                indelRef = md[x - 1][2]
                indelPos = md[x - 1][1] + 1
                while contLoop and md[x][0] is None:
                    indelSeq.append(md[x][2])
                    x += 1
                    if x == len(md):
                        contLoop = False
                indelDict.append(
                    f"{readChr}:"
                    f"{indelPos}:"
                    f"{''.join(indelSeq)}>{indelRef}"
                )
            else:
                x += 1
                if x == len(md):
                    contLoop = False
    return (indelDict)

def near_indel(inVar, inIndel):
    indelbins = inIndel.split(':')[2].split('>')
    indel_length = math.ceil(1.1 * max([len(indelbins[0]), len(indelbins[1])]))
    indelbins = inIndel.split(':')[:2]
    indel_start = int(indelbins[1]) - indel_length
    indel_stop = int(indelbins[1]) + indel_length
    varbins = inVar.split(':')[:2]
    if (
            varbins[0] == indelbins[0]
            and indel_start <= int(varbins[1]) <= indel_stop):
        return True
    else:
        return False

def is_indel(inVariant):
    varbins = inVariant.split(':')[2].split('>')
    if len(varbins[0]) != len(varbins[1]):
        return True
    else:
        return False

def extractVariant(inVarLine):
    if len(inVarLine.ref) == 1 and len(inVarLine.alts[0]) == 1:
        outStr = (
            f"{inVarLine.chrom}:{inVarLine.pos}"
            f":{refConvert[inVarLine.ref]}>"
            f"{compBase[inVarLine.alts[0]] if inVarLine.ref not in ('C', 'T') else inVarLine.alts[0]}"
        )
    else:
        outStr = (
            f"{inVarLine.chrom}:{inVarLine.pos}"
            f":{inVarLine.ref}>"
            f"{inVarLine.alts[0]}"
        )
    return outStr


def main():
    parser = ArgumentParser()
    parser.add_argument(
        '-i', '--inFile',
        action='store',
        dest='inFile',
        help='Input bam file.'
    )
    parser.add_argument(
        '-v', '--inVCF',
        action='store',
        dest='inVCF',
        help='Input VCF file for filtering.'
    )
    parser.add_argument(
        '-o', '--outPrefix',
        action='store',
        dest='outPrefix',
        help='Prefix for output files'
    )
    parser.add_argument(
        '-l', '--read_len',
        action='store',
        type=int,
        dest='rlen',
        help='The length of an individual read'
    )
    parser.add_argument(
        '--filter',
        action="append",
        dest="filters",
        default=[],
        help=(
            "Sets filtering behavior.  Can be provided multiple times to "
            "impose multiple filters.  'SNP' will filter out any "
            "variants marked as SNPs in the provided VCF file.  "
            "'INDEL' will filter out any variant located within indel_length "
            "bp of any indel in the provided VCF file.  'INCLUDE' "
            "will count only variants provided in the VCF file (with the "
            "exception of SNPs, if both 'INCLUDE' and 'SNP' are active).  "
            "All other arguments will be assumed to be filters in the 'FILTER'"
            " column of the VCF file, and will filter out variants matching "
            "that filter.  "
            ))
    parser.add_argument(
        '-g',
        action="store_true",
        dest="output_good",
        help="Output a file with reads below the threshold"
    )
    parser.add_argument(
        '-b',
        action="store_true",
        dest="output_bad",
        help="Output a file containing reads with less than mismatch count threshold."
    )
    parser.add_argument(
        '-t',
        action="store",
        dest="threshold",
        type=int,
        default=2,
        help="The threshold to determine if a read has too many mutations [2]."
    )
    parser.add_argument(
        '-c',
        action="store_true",
        dest="out_chart",
        help="Output a chart of the number of reads with various numbers of mismatches.  "
    )
    parser.add_argument(
        '--text_file',
        action="store_true",
        dest="output_text",
        help="Output a text file of the number of reads with various numbers of mismatches.  "
    )
    o = parser.parse_args()
    # Get list of SNPs from snp file
    inVCF = VariantFile(o.inVCF, 'r')
    filt_out = set()
    indels = set()
    variants = set()
    for line in inVCF:
        variant = extractVariant(line)
        if is_indel(variant) and "INDEL" in o.filters:
            indels.add(variant)
        else:
            filterVar = False
            for filt_iter in o.filters:
                if filt_iter not in ("INDEL", "INCLUDE"):
                    if filt_iter in line.filter:
                        filt_out.add(variant)
                        filterVar = True
            if "INCLUDE" in o.filters and not filterVar:
                variants.add(variant)

    myCounter = MismatchCounter(o.rlen, o.filters, filt_out, indels, variants)
    inBam = pysam.AlignmentFile(o.inFile, 'rb')

    # open output bam files
    if o.output_good:
        outGoodBam = pysam.AlignmentFile(
            f"{o.outPrefix}.goodReads.t{o.threshold}.bam",
            'wb',
            template=inBam
        )
    if o.output_bad:
        outBadBam = pysam.AlignmentFile(
            f"{o.outPrefix}.badReads.t{o.threshold}.bam",
            'wb',
            template=inBam
        )

    readCtr = 0
    numMismatchDict = defaultdict(int)
    for read in inBam:
        testSnps = myCounter.countRead(read)
        numMismatchDict[testSnps] += 1
        if testSnps > o.threshold:
            if o.output_bad:
                outBadBam.write(read)
        else:
            if o.output_good:
                outGoodBam.write(read)
        readCtr += 1
        if readCtr % 100000 == 0:
            print(f"{readCtr} reads processed...")
    # finish mutsPerRead output
    plotX = sorted(numMismatchDict.keys())
    plotY = [numMismatchDict[x] for x in plotX]
    if o.output_text:
        outTxt = open(f"{o.outPrefix}.mutsPerRead.txt", 'w')
        outTxt.write("Number_Mismatches\tNumber_Reads")
        for x in range(len(plotX)):
            outTxt.write(f"\n{plotX[x]}\t{plotY[x]}")
        outTxt.close()

    if o.out_chart:
        fig = plt.figure()
        plt.bar(x=plotX[o.threshold:], height=plotY[o.threshold:])
        plt.xlabel("Number of Mismatches")
        plt.ylabel("Number of Reads")
        plt.title(o.inFile)
        fig.savefig(f"{o.outPrefix}.mutsPerRead.png", bbox_inches='tight')
    if o.output_good:
        outGoodBam.close()
    if o.output_bad:
        outBadBam.close()

    myCounts = myCounter.getCounts()
    df = pd.DataFrame(myCounts)
    df['Base'] = range(1, len(df) + 1)
    df['Count_percent'] = df['Count'] / df['Count'].sum()
    df['Count_percent'] = df['Count_percent'] / df['Count_percent'].max()

    if o.out_chart:
        # Plot with Ns
        fig, ax = plt.subplots(figsize=(20, 14))
        margin_bottom = np.zeros(len(df['Base'].drop_duplicates()))
        colors = {"C>A": '#5abdeb', "C>G": '#050708', "C>T": '#d43c32', "T>A": '#cbcacb', "T>C": '#aacb72',
                  "T>G": '#e7c9c6', "C>N": "#4B0082", "T>N": "#4B0082", "": 'w'}

        for xIter in ["C>T", "C>A", "C>G", "T>A", "T>C", "T>G", "C>N", "T>N"]:
            df.plot.bar(x="Base", y=xIter, ax=ax, bottom=margin_bottom, color=colors[xIter]).set_ylabel("Count",
                                                                                                        fontsize=40)
            margin_bottom += df[xIter]

        labels = ax.get_xticklabels()

        for i, l in enumerate(labels):
            val = int(l.get_text())
            if val % 5 != 0:
                labels[i] = ''
            plt.gca().set_xticklabels(labels)

        ax.legend(fontsize=25)
        plt.xlabel("Cycle", fontsize=40)
        plt.yticks(fontsize=25)
        plt.xticks(fontsize=20)
        ax2 = ax.twinx()
        ax2.set_ylabel("Fraction of Total Reads", fontsize=40)
        ax2.set_title(o.inFile.split('/')[-1], fontsize=35)  # Not portable to Windows platform
        plt.plot(df["Base"], df['Count_percent'], color="black")
        plt.ylim(0, 1.02)
        plt.yticks(fontsize=25)
        fig.savefig(f"{o.outPrefix}_BasePerPosInclNs.png")

        # Plot without Ns
        fig, ax = plt.subplots(figsize=(20, 14))
        margin_bottom = np.zeros(len(df['Base'].drop_duplicates()))
        colors = {"C>A": '#5abdeb', "C>G": '#050708', "C>T": '#d43c32', "T>A": '#cbcacb', "T>C": '#aacb72',
                  "T>G": '#e7c9c6', "C>N": "#4B0082", "T>N": "#4B0082", "": 'w'}

        for xIter in ["C>T", "C>A", "C>G", "T>A", "T>C", "T>G"]:
            df.plot.bar(x="Base", y=xIter, ax=ax, bottom=margin_bottom, color=colors[xIter]).set_ylabel("Count",
                                                                                                        fontsize=40)
            margin_bottom += df[xIter]

        labels = ax.get_xticklabels()

        for i, l in enumerate(labels):
            val = int(l.get_text())
            if val % 5 != 0:
                labels[i] = ''
            plt.gca().set_xticklabels(labels)

        ax.legend(fontsize=25)
        plt.xlabel("Cycle", fontsize=40)
        plt.yticks(fontsize=25)
        plt.xticks(fontsize=20)
        ax2 = ax.twinx()
        ax2.set_ylabel("Fraction of Total Reads", fontsize=40)
        ax2.set_title(o.inFile.split('/')[-1], fontsize=35)  # Not portable to Windows platform
        plt.plot(df["Base"], df['Count_percent'], color="black")
        plt.ylim(0, 1.02)
        plt.yticks(fontsize=25)

        fig.savefig(f"{o.outPrefix}_BasePerPosWithoutNs.png")
    df = df[[
        "C>T", "C>A", "C>G",
        "T>A", "T>C", "T>G",
        "C>N", "T>N",
        "Count", "Base", "Count_percent"
    ]]
    df.to_csv(f"{o.outPrefix}_MutsPerCycle.dat.csv")


if __name__ == "__main__":
    main()
