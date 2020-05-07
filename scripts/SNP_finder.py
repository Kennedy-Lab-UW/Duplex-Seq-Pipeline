import logging
import argparse
from argparse import ArgumentParser
from VCF_Parser import *


def main():
    # Load Parameters
    parser = ArgumentParser()
    #  Input VCF file
    parser.add_argument(
        '-i', '--in_file',
        action='store',
        dest='Fin',
        help='An input single-sample VCF file to process',
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
    #  minimum depth for SNP detection
    parser.add_argument(
        '-d', '--min_depth',
        action='store',
        dest='min_depth',
        help='The minimum depth required to avoid being filtered',
        default=100,
        type=int
    )
    #  confidence
    parser.add_argument(
        '-c', '--confidence',
        action='store',
        dest='confidence',
        help='What level of confidence do you want to use in heterozygous SNP detection',
        default=0.05,
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

    # Load VCF file
    inVCF = VariantFile(o.Fin, 'r')
    # Check that this is a single-sample VCF file
    try:
        assert len(inVCF.samps) == 1
    except AssertionError:
        logging.error(f"Input VCF file {o.Fin} is not a single sample VCF.")
        raise
    sampName = inVCF.samps[0]

    # Add depth and snp filters to header
    inVCF.header.addLine(
        lineType="FILTER",
        label="low_depth",
        description=f"Sample depth at this locus is < {o.min_depth}."
    )
    inVCF.header.addLine(
        lineType="FILTER",
        label="SNP",
        description=f"This variant is probably a SNP, and not a somatic variant."
    )

    # Open SNP VCF
    outSNP = VariantFile(o.Fsnp, 'w', inVCF.header)

    # Open output VCF file
    outVCF = VariantFile(o.Fout, 'w', inVCF.header)

    # For variant in VCF file
    for varLine in inVCF:
        #  If depth < minDepth
        if int(varLine.samples[sampName]["DP"]) < o.min_depth:
            #   Add low_depth filter
            varLine.add_filter("low_depth")
        if float(varLine.samples[sampName]["AF"].split(',')[1]) >= 0.4:
            # Mark as SNP
            varLine.add_filter("SNP")
            # Write to SNP file
            outSNP.writeline(varLine)
        #  write mutation to output file
        outVCF.writeline(varLine)
    # Close files
    inVCF.close()
    outSNP.close()
    outVCF.close()
    # DONE!


if __name__ == "__main__":
    main()
