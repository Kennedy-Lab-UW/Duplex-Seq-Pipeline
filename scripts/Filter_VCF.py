import logging
import argparse
from argparse import ArgumentParser
import math
from collections import namedtuple
import pysam
import regex
from VCF_Parser import VariantFile
from BedParser import Bed_File

def extractVariant(inVarLine):
    outStr = (
        f"{inVarLine.chrom}:{inVarLine.pos}"
        f":{inVarLine.ref}>"
        f"{inVarLine.alts[0]}"
    )
    trueLen = max([len(inVarLine.ref), len(inVarLine.alts[0])])
    truePos = inVarLine.pos
    if len(inVarLine.ref) > 1 and len(inVarLine.alts[0]) == 1:
        varType = "del"
    elif len(inVarLine.ref) == 1 and len(inVarLine.alts[0]) > 1:
        varType = "ins"
    else:
        varType = "complex"
    #varType = "ins" if len(inVarLine.ref) > 1 else "del" if 
    lenOfVar = math.ceil(1.1 * trueLen)
    outStart = inVarLine.pos - lenOfVar
    outStop = inVarLine.pos + lenOfVar

    return outStr, inVarLine.chrom, outStart, outStop, truePos, trueLen, varType

indelRegion = namedtuple(
    "IndelRegion",
    ["chrom",
     "start",
     "stop",
     "truePos",
     "trueLen",
     "varType"])

def main():
    """main program"""
    # Arguments
    parser = ArgumentParser()
    parser.add_argument(
        '-i', '--in_vcf',
        action='store',
        dest='in_vcf',
        help='An input VCF file to mask.',
        required=True)
    parser.add_argument(
        '-o', '--out_vcf',
        action='store',
        dest='out_vcf',
        help='A name for the output (masked) VCF file',
        required=True)
    parser.add_argument(
        '-s', '--snp_file',
        action='store',
        dest='Fsnp',
        help='An name for the output VCF file including only SNPs.',
        required=True)
    parser.add_argument(
        '-b', '--mask_bed',
        action='store',
        dest='bed',
        default=None,
        help='A bed file with the regions to be masked.')
    parser.add_argument(
        '-d', '--min_depth',
        action='store',
        dest='min_depth',
        help='The minimum depth required to avoid being filtered.',
        default=100,
        type=int)
    parser.add_argument(
        '--n_lim', 
        action='store',
        dest='n_lim',
        type=float,
        help='The maximum proportion of Ns to allow to avoid being filtered.',
        default=0.1)
    parser.add_argument(
        '-t', '--snp_threshold',
        action='store',
        dest='snp_threshold',
        help=(f'The threshold for a variant to be marked as a SNP.  '
              f'Any mutation with a MAF >= this will be marked as a SNP.'),
        default=0.4,
        type=float)
    parser.add_argument(
        '--cluster_dist', 
        action='store', 
        dest='cluster_dist',
        type=int,
        default=10,
        help='How far from a variant to check for clustering.')
    parser.add_argument(
        "--max_indel_size", 
        action="store",
        dest="max_indel_size",
        type=int,
        default=20,
        help="The maximum size of indel to keep without filtering")    
    # logLevel (hidden argument)
    parser.add_argument(
        '--logLevel',
        action="store",
        dest="logLvl",
        default="INFO",
        help=argparse.SUPPRESS,
        choices={'INFO', 'DEBUG', 'WARNING', 'ERROR', 'CRITICAL'})
    o = parser.parse_args()

    # Setup logging
    numeric_level = getattr(logging, o.logLvl.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % o.logLvl.upper())
    logging.basicConfig(
        format='%(levelname)s: %(message)s',
        level=numeric_level)

    logging.info(f"Masking VCF file {o.in_vcf}")
    logging.info(f"Using bed file {o.bed}")
    logging.info(f"Writing output to VCF file {o.out_vcf}")
    # Read input VCF
    in_vcf = VariantFile(o.in_vcf, 'r')
    # Extract header
    my_header = in_vcf.header
    # Get sample name
    o.sampName = in_vcf.samps[0]
    # Read variants
    myVariants = [x for x in in_vcf]
    mySnps = []
    indels = {}

    # add SNP filter to header
    my_header.addLine(
        lineType="FILTER",
        label="SNP",
        description=f"This variant is probably a SNP (MAF >= {o.snp_threshold}), and not a somatic variant.")
    # add low-depth filter to header
    my_header.addLine(
        lineType="FILTER",
        label="low_depth",
        description=f"Sample depth at this locus is < {o.min_depth}.")
    # add near-indel filter to header
    my_header.addLine(
        lineType="FILTER",
        label="near_indel",
        description=f"This variant may be suspecious because it is near an indel.")
    my_header.addLine(
        lineType="FILTER",
        label="near_complex",
        description=f"This variant may be suspecious because it is near a complex variant (MNP, or insertion and deletion in the same area, or indel and SNP in the same area).")
    # add cluster filter to header
    my_header.addLine(
        lineType="FILTER",
        label="clustered",
        description=f"This variant is part of a cluster of variants (nearest variant is within {o.cluster_dist} bp).")
    my_header.addLine(
        lineType="FILTER",
        label="high_Ns",
        description=f"This variant has a high proportion of Ns (> {o.n_lim}).")
    my_header.addLine(
        lineType="FILTER",
        label="long_indel",
        description=f"This variant is an improbably long indel (len(indel) > {o.max_indel_size}).")
    # Add "ins", "del", and "complex" filters (mostly just markers that the variants have these types)
    my_header.addLine(
        lineType="FILTER",
        label="ins",
        description=f"This variant is an insertion.")
    my_header.addLine(
        lineType="FILTER",
        label="del",
        description=f"This variant is a deletion.")
    my_header.addLine(
        lineType="FILTER",
        label="complex",
        description=f"This variant is a complex variant.")
    # Add "masked" filter to header
    my_header.addLine(
        "FILTER",
        "masked",
        description=f"This position has been masked due to overlap with {o.bed}.")
    # Add "complex_near_SNP" filter to header
    my_header.addLine(
        "FILTER",
        "complex_near_SNP",
        description="This variant is a complex variant that overlaps a SNP"
    )
    # Get set of SNP loci and indel loci, and add low_depth and SNP filters
    for varLine in myVariants:
        if int(varLine.samples[o.sampName]["DP"]) < o.min_depth:
            #   Add low_depth filter
            varLine.add_filter("low_depth")
        if float(varLine.samples[o.sampName]["AF"].split(',')[1]) >= o.snp_threshold:
            # Mark as SNP
            varLine.add_filter("SNP")
        if (float(varLine.samples[o.sampName]["NC"]) / 
                (float(varLine.samples[o.sampName]["DP"]) + 
                 float(varLine.samples[o.sampName]["NC"]) ) > o.n_lim):
            varLine.add_filter("high_Ns")
        if len(varLine.ref) > 1 or max([len(x) for x in varLine.alts]) > 1:
            myVar = extractVariant(varLine)
            indels[myVar[0]] = indelRegion(myVar[1], myVar[2], myVar[3],
                                           myVar[4], myVar[5], myVar[6])
            varLine.add_filter(myVar[6])
            if myVar[5] > o.max_indel_size and myVar[6] is not "complex":
                varLine.add_filter("long_indel")
    myVars2 = []
    for vcfIter in range(len(myVariants)):
        nearby_variants = []
        is_ins = False
        is_del = False
        is_complex = False
        isIndel = False
        if  "long_indel" not in myVariants[vcfIter].filter:
            if "ins" in myVariants[vcfIter].filter:
                is_ins = True
            if "del" in myVariants[vcfIter].filter:
                is_del = True
            if "complex" in myVariants[vcfIter].filter:
                is_complex = True            
        if is_ins or is_del or is_complex:
            thisVariant = indelRegion(*extractVariant(myVariants[vcfIter])[1:])
            isIndel = True
        for vcfIter2 in range(len(myVariants)):
            if vcfIter == vcfIter2:
                continue
            nearIndel = False
            nearComplex = False
            compNearSnp = False
            # check for overlapping complex variants near a SNP
            if ("complex" in myVariants[vcfIter2].filter and 
                    "SNP" in myVariants[vcfIter].filter and
                    myVariants[vcfIter].chrom == myVariants[vcfIter2].chrom and
                    myVariants[vcfIter].pos >= myVariants[vcfIter2].pos and 
                    myVariants[vcfIter].pos <= myVariants[vcfIter2].pos + len(myVariants[vcfIter2].ref)):
                compNearSnp = True
            # check for variant near an indel
            if isIndel:
                if (thisVariant.chrom == myVariants[vcfIter2].chrom 
                        and thisVariant.start <= myVariants[vcfIter2].pos <= thisVariant.stop
                        and not is_complex):
                    nearIndel = True
                if (thisVariant.chrom == myVariants[vcfIter2].chrom 
                        and thisVariant.truePos <= myVariants[vcfIter2].pos
                        and myVariants[vcfIter2].pos <= thisVariant.truePos + thisVariant.trueLen
                        and is_complex):
                    nearComplex = True
            if compNearSnp and "SNP" not in myVariants[vcfIter].filter: 
                myVariants[vcfIter2].add_filter("complex_near_SNP")
            if nearIndel and "SNP" not in myVariants[vcfIter2].filter:
                myVariants[vcfIter2].add_filter("near_indel")
            if nearComplex and "SNP" not in myVariants[vcfIter2].filter:
                myVariants[vcfIter2].add_filter("near_complex")
            
            # add clustered variant filter
            if  ("SNP" not in myVariants[vcfIter].filter
                    and myVariants[vcfIter].chrom == myVariants[vcfIter2].chrom
                    and abs(myVariants[vcfIter].pos - myVariants[vcfIter2].pos) <= o.cluster_dist
                    and "SNP" not in myVariants[vcfIter2].filter):
                nearby_variants.append(vcfIter2)
        if len(nearby_variants) > 0:
            myVariants[vcfIter].add_filter("clustered")
      
    # Open output VCF file
    out_vcf = VariantFile(o.out_vcf, 'w', my_header)
    # Open bed file, and extract regions
    if o.bed is not None:
        mask_bed = [x for x in Bed_File(o.bed)]
    else:
        mask_bed = []
    # Activate counters
    lines_read = 0
    lines_written = 0
    lines_masked = 0

    # Iterate over VCF lines
    for line in myVariants:
        lines_read += 1
        # Test if this line overlaps a region in the masking bed file
        mask_line = False
        for region in mask_bed:
            if region.contains(line.chrom, line.pos - 1):
                mask_line = True
                # stop iteration over regions if line is masked
                break
        # Apply masking
        if mask_line:
            lines_masked += 1
            line.add_filter("masked")
        # Check if line is a SNP
        if "SNP" in line.filter:
            mySnps.append(line)
        # Write line
        out_vcf.writeline(line)
        lines_written += 1
    #write output  SNPs file
    outSNPS = VariantFile(o.Fsnp, 'w', my_header)
    for vcfLine in mySnps:
        outSNPS.writeline(vcfLine)
    outSNPS.close()

    # close vcf flies
    in_vcf.close()
    out_vcf.close()
    logging.info(f"Read {lines_read} VCF lines")
    logging.info(f"Masked {lines_masked} lines")
    logging.info(f"Wrote {lines_written} lines")
    logging.info("DONE")

if __name__ == "__main__":
    main()
