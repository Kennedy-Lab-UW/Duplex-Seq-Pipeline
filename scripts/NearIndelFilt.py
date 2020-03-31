import argparse
from argparse import ArgumentParser
import collections
from collections import defaultdict, namedtuple
import pysam
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from VCF_Parser import *
import sys
import math


def extractVariant(inVarLine):
    outStr = (
        f"{inVarLine.chrom}:{inVarLine.pos}"
        f":{inVarLine.ref}>"
        f"{inVarLine.alts[0]}"
        )
    lenOfVar = math.ceil(1.1 * max([len(inVarLine.ref),len(inVarLine.alts[0])]))
    outStart = inVarLine.pos - lenOfVar
    outStop  = inVarLine.pos + lenOfVar
    return(outStr,inVarLine.chrom, outStart,outStop)
    
indelRegion = namedtuple(
    "IndelRegion", 
    [
        "chrom",
        "start",
        "stop"
        ]
    )
        
def main():
    parser = ArgumentParser()
    parser.add_argument(
        '-i','--inFile', 
        action='store', 
        dest='inFile', 
        help='Input VCF file.'
        )
    parser.add_argument(
        '-v', '--inIndels', 
        action='store', 
        dest='inIndels', 
        help='Input Indel VCF file.'
        )
    parser.add_argument(
        '-o', '--outFile', 
        action='store', 
        dest='outFile', 
        help='Name for output VCF file'
        )
    o = parser.parse_args()
    # Get list of SNPs from snp file
    inIndels = VariantFile(o.inIndels, 'r')
    indels = {}
    for line in inIndels:
        if "SNP" in line.filter:
            myVar = extractVariant(line)
            indels[myVar[0]] = indelRegion(myVar[1], myVar[2],myVar[3])
    inIndels.close()
    # ~ print(indels)
    inVCF = VariantFile(o.inFile,'r')
    inVCF.header.addLine(
        lineType="FILTER", 
        label="near_indel", 
        description=f"This variant may be suspecious because it is near an indel."
        )
    outVCF = VariantFile(o.outFile, 'w', inVCF.header)
    recordsRead = 0
    recordsWritten = 0
    for vcfLine in inVCF:
        recordsRead += 1
        nearIndel = False
        for indelIter in indels:
            if (
                    vcfLine.chrom == indels[indelIter].chrom
                    and vcfLine.pos >= indels[indelIter].start
                    and vcfLine.pos <= indels[indelIter].stop 
                    and len(vcfLine.ref) == 1
                    and len(vcfLine.alts[0]) == 1
                    ):
                nearIndel = True
            if nearIndel:
                vcfLine.add_filter("near_indel")
        outVCF.writeline(vcfLine)
        recordsWritten += 1
    # ~ print(f"{recordsRead} read\n{recordsWritten} written")
    inVCF.close()
    outVCF.close()

if __name__ == "__main__":
    main()
