import argparse
from argparse import ArgumentParser
import logging
import sys
from collections import defaultdict

import pysam
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from VCF_Parser import *

refConvert = {"C": "C", "G": "C", "T": "T", "A": "T"}
compBase = {"C": "G", "G": "C", "T": "A", "A": "T", "N": "N"}

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

def getNonOverlap(inRead, use_R2=False):
    inRead.query_length
    inRead.is_read1
    if (2 * inRead.query_length >= inRead.template_length 
            and inRead.template_length > 0
            and not inRead.is_reverse):
        if ((inRead.is_read1 and not use_R2)
                or (inRead.is_read2 and use_R2)):
            return {"start": 0, "end": inRead.query_length}
        elif ((inRead.is_read2 and not use_R2)
                or (inRead.is_read1 and use_R2)):
            return {"start": 0, "end": inRead.query_length*2-abs(inRead.template_length)}
    elif (2 * inRead.query_length >= -inRead.template_length 
            and inRead.template_length < 0
            and inRead.is_reverse):
        if ((inRead.is_read1 and not use_R2)
                or (inRead.is_read2 and use_R2)):
            return {"start": 0, "end": inRead.query_length}
        elif ((inRead.is_read2 and not use_R2)
                or (inRead.is_read1 and use_R2)):
            return {"start": 0, "end": inRead.query_length*2-abs(inRead.template_length)}
    else:
        return {"start": 0, "end": inRead.query_length}

def main():
    # Load parameters
    parser = ArgumentParser()
    parser.add_argument(
        '--inVCF',
        action='store',
        dest='inVCF',
        help='Input VCF file.',
        required=True
    )
    parser.add_argument(
        '--inBam',
        action='store',
        dest='inBAM',
        help='Input BAM file.',
        required=True
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
        '--clonalityCutoff',
        action='store',
        type=float,
        dest='clonal',
        help='The maximum clonality to use for muts_by_read_position.',
        default=.1
    )
    parser.add_argument(
        '--depthCutoff',
        action='store',
        type=int,
        dest='depth',
        help='The minimum depth to use for muts_by_read_position.',
        default=100
    )
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
    parser.add_argument(
        '-5',
        action="store",
        type=int,
        dest="clip5",
        help="Amount of 5' clipping"
    )
    parser.add_argument(
        '-3',
        action="store",
        type=int,
        dest="clip3",
        help="Amount of 3' clipping"
    )
    parser.add_argument(
        '--useR2',
        action="store_true",
        dest="useR2",
        help="If the two reads overlap, use R2 over R1."
    )
    parser.add_argument(
        '--unique',
        action="store_true",
        dest="unique",
        help="Only count each mutation once."
    )
    parser.add_argument(
        '--logLevel',
        action="store",
        dest="logLvl",
        default="INFO",
        help=argparse.SUPPRESS,
        choices={'INFO', 'DEBUG', 'WARNING', 'ERROR', 'CRITICAL'}
    )
    o = parser.parse_args()
    
    # Setup logging
    numeric_level = getattr(logging, o.logLvl.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % o.logLvl.upper())
    logging.basicConfig(
        format='%(levelname)s: %(message)s',
        level=numeric_level,
    )
    
    # Load variant position data
    # Establish data storage
    mismatch_counts = [
        {"Base": i+1, "C>T": 0, "C>A": 0, "C>G": 0, "T>A": 0, "T>C": 0, "T>G": 0, "C>N": 0, "T>N": 0, "Count": 0}
        for i in range(o.rlen)]
    inVCF = VariantFile(o.inVCF, 'r')
    sampID = inVCF.samps[0]
    for varLine in inVCF:
        # check the clonality of the variant
        varClonal = float(varLine.samples[sampID]['AF'].split(',')[1])
        varAD = int(varLine.samples[sampID]['AD'].split(',')[1])
        varDP = int(varLine.samples[sampID]['DP'])
        varPM = int(float(varLine.samples[sampID]['PM']))
        if (((0 < varClonal <= o.clonal or varAD == 1) and varDP >= o.depth)
                and len(varLine.ref) == 1 and len(varLine.alts[0]) == 1):
            changeID = (
                f"{refConvert[varLine.ref]}>"
                f"{compBase[varLine.alts[0]] if varLine.ref not in ('C', 'T') else varLine.alts[0]}"
            )
            
            mismatch_counts[varPM - 1][changeID] += 1 if o.unique else varAD
    
    # Load read data
    inBam = pysam.AlignmentFile(o.inBAM, 'rb')
    readCounter = 0
    for read in inBam:
        # apply clipping
        bounds = {"start": o.clip5, "end": read.query_length - o.clip3}
        
        # check for overlap
        overlap = getNonOverlap(read, o.useR2)
        # merge clipping
        finalBounds = {"start":max(bounds["start"], overlap["start"]), 
            "end": min(bounds["end"], overlap["end"])}
        for xIter in range(finalBounds["start"],finalBounds["end"]):
            mismatch_counts[xIter]["Count"] += 1
        readCounter += 1
        if readCounter % 100000 == 0:
            logging.info(f"{readCounter} reads processed...")
    # Write outputs:
    df = pd.DataFrame(mismatch_counts)
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
        ax2.set_title(o.inVCF.split('/')[-1], fontsize=35)  # Not portable to Windows platform
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
        ax2.set_title(o.inVCF.split('/')[-1], fontsize=35)  # Not portable to Windows platform
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