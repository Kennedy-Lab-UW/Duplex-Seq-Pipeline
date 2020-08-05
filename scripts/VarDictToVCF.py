import sys
import logging
import argparse
from argparse import ArgumentParser
import math
from collections import namedtuple
from VCF_Parser import *
import pandas as pd
import pysam
import regex


def varDictLine2vcfLine(inVarDictLine, inSampName):
    recordStr = (
        f"{inVarDictLine['Chr']}\t"
        f"{inVarDictLine['Start']}\t"
        f".\t"  # id
        f"{inVarDictLine['Ref']}\t"
        f"{inVarDictLine['Alt']}\t"
        f".\t"  # qual
        f".\t"  # filter
        f"GENE={inVarDictLine['Gene']};"
        f"MSI={inVarDictLine['MSI']};"
        f"MSINT={inVarDictLine['MSI_NT']};"
        f"5pFlank={inVarDictLine['5pFlankSeq']};"
        f"3pFlank={inVarDictLine['3pFlankSeq']};"
        f"VARTYPE={inVarDictLine['VarType']}\t"  # info
        f"AD:DP:AF:NC:FC:RC:PM:QM:MQ:NM\t"  # format
        f"{inVarDictLine['RefCount']},{inVarDictLine['AltDepth']}:"
        f"{inVarDictLine['Depth']}:"
        f"{inVarDictLine['RefFract']:.2E},{inVarDictLine['AltFract']:.2E}:"
        f"{inVarDictLine['NC']}:"
        f"{inVarDictLine['RefFwdReads']},{inVarDictLine['AltFwdReads']}:"
        f"{inVarDictLine['RefRevReads']},{inVarDictLine['AltRevReads']}:"
        f"{inVarDictLine['PMean']}:"
        f"{inVarDictLine['QMean']}:"
        f"{inVarDictLine['MQ']}:"
        f"{inVarDictLine['NM']}"
    )
    return VariantRecord(recordStr, [inSampName])


def extractVariant(inVarLine):
    outStr = (
        f"{inVarLine.chrom}:{inVarLine.pos}"
        f":{inVarLine.ref}>"
        f"{inVarLine.alts[0]}"
    )
    lenOfVar = math.ceil(1.1 * max([len(inVarLine.ref), len(inVarLine.alts[0])]))
    outStart = inVarLine.pos - lenOfVar
    outStop = inVarLine.pos + lenOfVar
    return outStr, inVarLine.chrom, outStart, outStop


indelRegion = namedtuple(
    "IndelRegion",
    ["chrom",
     "start",
     "stop"])


def main():
    parser = ArgumentParser()
    #  Input VarDict file
    parser.add_argument(
        '-i', '--in_file',
        action='store',
        dest='Fin',
        help=(f"An input single-sample VarDict output file without Ns "
              f"included in the depth to process"),
        required=True
    )
    parser.add_argument(
        '-n', '--in_Nfile',
        action='store',
        dest='Nin',
        help=(f"An input single-sample VarDict output file with Ns "
              f"included in the depth to process"),
        required=True
    )
    parser.add_argument(
        '-b', '--in_bam_file',
        action='store',
        dest='Bin',
        help=f"An input bam file used to generate both VarDict inputs",
        required=True
    )
    #  Output VCF file
    parser.add_argument(
        '-o', '--out_file',
        action='store',
        dest='Fout',
        help='An name for the output VCF file including depth and SNP filtering',
        required=True
    )
    #  Output SNP file
    parser.add_argument(
        '-s', '--snp_file',
        action='store',
        dest='Fsnp',
        help='An name for the output VCF file including only SNPs',
        required=True
    )
    #  minimum depth for good variant detection
    parser.add_argument(
        '-d', '--min_depth',
        action='store',
        dest='min_depth',
        help='The minimum depth required to avoid being filtered',
        default=100,
        type=int
    )
    parser.add_argument(
        '--samp_name',
        action='store',
        dest='sampName',
        help='A name for this sample.  ',
        default="SAMPLE"
    )
    parser.add_argument(
        '-t', '--snp_threshold',
        action='store',
        dest='snp_threshold',
        help=(f'The threshold for a variant to be marked as a SNP.  '
              f'Any mutation with a MAF >= this will be marked as a SNP.  '),
        default=0.4,
        type=float
    )
    # logLevel (hidden argument)
    parser.add_argument(
        '--logLevel',
        action="store",
        dest="logLvl",
        default="INFO",
        help=argparse.SUPPRESS,
        choices={'INFO', 'DEBUG', 'WARNING', 'ERROR', 'CRITICAL'}
    )
    # parse arguments
    o = parser.parse_args()

    # Setup logging
    numeric_level = getattr(logging, o.logLvl.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % o.logLvl.upper())
    logging.basicConfig(
        format='%(levelname)s: %(message)s',
        level=numeric_level,
    )
    # SWITCH THIS LINE IN BEFORE DEPLOYMENT; it will get the current version from VERSION
    # pipelineVersion=open(f"{sys.path[0]}/../VERSION", 'r').readline().strip()
    pipelineVersion = "2.0.0b"
    cmd = " ".join(sys.argv)
    # create VCF header
    myHeader = SamHeaderToVcfHeader(
        pysam.AlignmentFile(o.Bin, "rb").header,
        o.sampName,
        "VarDictToVCF",
        pipelineVersion,
        cmd)
    # Numbers: #, A, R, G, .
    # Types (Format): Integer, Float, Character, String
    # Types (Info): Integer, Float, Flag, Character, String
    vardictVersion = "1.7.0"
    myHeader.addLine(lineType="INFO",
                     label="GENE",
                     number="1",
                     Type="Character",
                     description="Target name from the provided bed file",
                     source="VarDict-Java",
                     vNum=vardictVersion)
    myHeader.addLine(lineType="INFO",
                     label="MSI",
                     number="1",
                     Type="Integer",
                     description="MicroSatellite. > 1 indicates MSI",
                     source="VarDict-Java",
                     vNum=vardictVersion)
    myHeader.addLine(lineType="INFO",
                     label="MSINT",
                     number="1",
                     Type="Integer",
                     description="MicroSatellite unit length in bp`",
                     source="VarDict-Java",
                     vNum=vardictVersion)
    myHeader.addLine(lineType="INFO",
                     label="5pFlank",
                     number="1",
                     Type="Character",
                     description="Neighboring reference sequence to 5\' end",
                     source="VarDict-Java",
                     vNum=vardictVersion)
    myHeader.addLine(lineType="INFO",
                     label="3pFlank",
                     number="1",
                     Type="Character",
                     description="Neighboring reference sequence to 3\' end",
                     source="VarDict-Java",
                     vNum=vardictVersion)
    myHeader.addLine(lineType="INFO",
                     label="VARTYPE",
                     number="1",
                     Type="Character",
                     description="Variant type",
                     source="VarDict-Java",
                     vNum=vardictVersion)
    myHeader.addLine(lineType="INFO",
                     label="VarDesc",
                     number="1",
                     Type="Character",
                     description="Variant description string, from VarDict Genotype field",
                     source="VarDict-Java",
                     vNum=vardictVersion)
    myHeader.addLine(lineType="FORMAT",
                     label="AD",
                     number="R",
                     Type="Integer",
                     description="Depth contribution of each allele")
    myHeader.addLine(lineType="FORMAT",
                     label="DP",
                     number="1",
                     Type="Integer",
                     description="Total depth")
    myHeader.addLine(lineType="FORMAT",
                     label="AF",
                     number="R",
                     Type="Float",
                     description="Fraction of total depth accounted for by each allele")
    myHeader.addLine(lineType="FORMAT",
                     label="NC",
                     number="1",
                     Type="Integer",
                     description="Number of Ns at this position")
    myHeader.addLine(lineType="FORMAT",
                     label="AFC",
                     number="A",
                     Type="Integer",
                     description="The number of reads for each allele that are forward mapping")
    myHeader.addLine(lineType="FORMAT",
                     label="ARC",
                     number="A",
                     Type="Integer",
                     description="The number of reads for each allele that are reverse mapping")
    myHeader.addLine(lineType="FORMAT",
                     label="PM",
                     number="1",
                     Type="Float",
                     description="Mean position of the alt allele in the read")
    myHeader.addLine(lineType="FORMAT",
                     label="QM",
                     number="1",
                     Type="Float",
                     description="Mean quality score of data at this position")
    myHeader.addLine(lineType="FORMAT",
                     label="MQ",
                     number="1",
                     Type="Float",
                     description="Mean mapping quality of reads covering this position")
    myHeader.addLine(lineType="FORMAT",
                     label="NM",
                     number="1",
                     Type="Float",
                     description="Average number of mismatches for reads containing the ALT allele")
    # read in non-Ns vardict file and sort it
    try:
        vardict_vars = pd.read_csv(o.Fin, sep='\t', header=0).sort_values(['Chr', 'Start'])
    except KeyError:
        message = f"in_file {o.Fin} appears to lack a header row.  Did you generate this file using -h?"
        logging.error(message)
        raise
    vardict_vars["varID"] = vardict_vars["Chr"] + ":" + vardict_vars["Start"].map(str) + "-" + vardict_vars["End"].map(
        str)
    # read in Ns vardict file and sort it
    try:
        vardict_Nvars = pd.read_csv(o.Nin, sep='\t', header=0).sort_values(['Chr', 'Start'])
    except KeyError:
        message = f"in_Nfile {o.Nin} appears to lack a header row.  Did you generate this file using -h?"
        logging.error(message)
        raise
    vardict_Nvars["varID"] = vardict_Nvars["Chr"] + ":" + vardict_Nvars["Start"].map(str) + "-" + vardict_Nvars[
        "End"].map(str)
    # check that the two vardict files contain the same positions
    assert all(vardict_vars["varID"] == vardict_Nvars["varID"])

    # make Ns column in non-Ns vardict data
    vardict_vars["NC"] = vardict_Nvars["Depth"] - vardict_vars["Depth"]

    # calculate AF columns in vardict data
    vardict_vars["RefCount"] = vardict_vars["Depth"] - vardict_vars["AltDepth"]
    vardict_vars["RefFract"] = vardict_vars["RefCount"] / vardict_vars["Depth"]
    vardict_vars["AltFract"] = vardict_vars["AltDepth"] / vardict_vars["Depth"]

    # create array of VCF lines for vardict lines
    myVariants = []
    for index, rowIter in vardict_vars.iterrows():
        if (len(regex.findall(r'^[ACGTN]+$',str(rowIter['Ref']))) == 1 
                and len(regex.findall(r'^[ACGTN]+$',str(rowIter['Alt']))) == 1):
            if (rowIter["AltDepth"] > 0
                    and 'N' not in rowIter['Ref']
                    and 'N' not in rowIter['Alt']):
                myVariants.append(varDictLine2vcfLine(rowIter, o.sampName))
    mySnps = []
    indels = {}
    # add SNP filter to header
    myHeader.addLine(
        lineType="FILTER",
        label="SNP",
        description=f"This variant is probably a SNP (MAF >= {o.snp_threshold}), and not a somatic variant.")
    # add low-depth filter to header
    myHeader.addLine(
        lineType="FILTER",
        label="low_depth",
        description=f"Sample depth at this locus is < {o.min_depth}.")
    # add near-indel filter to header
    myHeader.addLine(
        lineType="FILTER",
        label="near_indel",
        description=f"This variant may be suspecious because it is near an indel.")
    # add cluster filter to header
    myHeader.addLine(
        lineType="FILTER",
        label="clustered",
        description=f"This variant is part of a cluster of variants (nearest variant is within 10 bp).")
    # Get set of SNP loci and indel loci, and add low_depth and SNP filters
    for varLine in myVariants:
        if int(varLine.samples[o.sampName]["DP"]) < o.min_depth:
            #   Add low_depth filter
            varLine.add_filter("low_depth")
        if float(varLine.samples[o.sampName]["AF"].split(',')[1]) >= o.snp_threshold:
            # Mark as SNP
            varLine.add_filter("SNP")
            # Write to SNP file
            mySnps.append(varLine)
        if len(varLine.ref) > 1 or max([len(x) for x in varLine.alts]) > 1:
            myVar = extractVariant(varLine)
            indels[myVar[0]] = indelRegion(myVar[1], myVar[2], myVar[3])

    # add near_indel filter
    for vcfLine in myVariants:
        nearIndel = False
        for indelIter in indels:
            if (
                    vcfLine.chrom == indels[indelIter].chrom
                    and indels[indelIter].start <= vcfLine.pos <= indels[indelIter].stop
                    and len(vcfLine.ref) == 1
                    and len(vcfLine.alts[0]) == 1):
                nearIndel = True
            if nearIndel:
                vcfLine.add_filter("near_indel")
    # add clustered variant filter
    for vcfIter in range(len(myVariants)):
        clusteredVariant = False
        if vcfIter == 0 and len(myVariants) > 1:
            if (
                    myVariants[vcfIter].chrom == myVariants[vcfIter + 1].chrom
                    and myVariants[vcfIter].pos + 10 >= myVariants[vcfIter + 1].pos):
                clusteredVariant = True
        elif vcfIter == len(myVariants) - 1 and len(myVariants) > 1:
            if (
                    myVariants[vcfIter].chrom == myVariants[vcfIter - 1].chrom
                    and myVariants[vcfIter].pos - 10 <= myVariants[vcfIter - 1].pos):
                clusteredVariant = True
        else:
            if ((
                    myVariants[vcfIter].chrom == myVariants[vcfIter - 1].chrom
                    and myVariants[vcfIter].pos - 10 <= myVariants[vcfIter - 1].pos) or (
                    myVariants[vcfIter].chrom == myVariants[vcfIter + 1].chrom
                    and myVariants[vcfIter].pos + 10 >= myVariants[vcfIter + 1].pos)):
                clusteredVariant = True
        if clusteredVariant:
            myVariants[vcfIter].add_filter("clustered")
    # Create output file
    outVCF = VariantFile(o.Fout, 'w', myHeader)
    # write lines
    for vcfLine in myVariants:
        outVCF.writeline(vcfLine)
    outVCF.close()
    outSNPS = VariantFile(o.Fsnp, 'w', myHeader)
    for vcfLine in mySnps:
        outSNPS.writeline(vcfLine)
    outSNPS.close()


if __name__ == "__main__":
    main()
