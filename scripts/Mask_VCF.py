import logging
import argparse
from argparse import ArgumentParser
from VCF_Parser import VariantFile
from BedParser import Bed_File

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
        '-b', '--mask_bed',
        action='store',
        dest='bed',
        help='A bed file with the regions to be masked.',
        required=True)
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
    # Add "masked" filter to header
    my_header = in_vcf.header
    my_header.addLine(
        "FILTER",
        "masked",
        description=f"This position has been masked due to overlap with {o.bed}.")
    # Open output VCF file
    out_vcf = VariantFile(o.out_vcf, 'w', my_header)
    # Open bed file, and extract regions
    mask_bed = [x for x in Bed_File(o.bed)]
    # Activate counters
    lines_read = 0
    lines_written = 0
    lines_masked = 0

    # Iterate over VCF lines
    for line in in_vcf:
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
        # Write line
        out_vcf.writeline(line)
        lines_written += 1
    # close vcf flies
    in_vcf.close()
    out_vcf.close()
    logging.info(f"Read {lines_read} VCF lines")
    logging.info(f"Masked {lines_masked} lines")
    logging.info(f"Wrote {lines_written} lines")
    logging.info("DONE")

if __name__ == "__main__":
    main()
